/*
 * Simple Dirichlet/Pitman-Yor Process sampling test
 * Copyright (C) 2011-2012 Wray Buntine and Lan Du
 *           (C) 2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *         Lan Du (lan.du@nicta.com.au)
 *     
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yaps.h"
#include "stable.h"
#include "lgamma.h"
#include "digamma.h"
#include "psample.h"
/*
 *    the RNG is GSL, so allocation, reallocation and freeing
 *    all done with GSL calls, but all sampling using rng_*()
 *    interface here
 */
#include "srng.h"

/*
 *   Two options for sampling hyperparameters:
 *          Adaptive rejection sampling (ARS)
 *              ->  complex, we use Wally Gilks' C library
 *          Slice sampling
 *              ->  very simple, but needs to be warmed up
 *          For sampling b,a, these selected with a switch in psample.h
 */
#include "arms.h"


/*==================================================
 * global variables
 *
 *     This is a very primitive setup to test a single
 *     Pitman-Yor.
 *==================================================*/

rngp_t rng = NULL;
int verbose = 0;

/*
 *    dimensions and constants
 */
//    H() is discrete distribution with this dimension
#define DIM 20
//    number of multinomials sampled
#define NUMMN 3
//    array limits
#define MAXSTAB 1000
#define MAXTAB 1000
#define MAXDATA 10000

/*
 *   parameters and hyper parameters
 */
float  apar = 0.0, bpar = 1.0, binit;
#define PB_shape 1.1
#define PB_scale DIM
//    base distribution
float H[DIM];

/*
 *   data
 */
int data[MAXDATA];
scnt_int N[NUMMN]; 
//    statistics
scnt_int *n[NUMMN];
scnt_int n_data[DIM*NUMMN];
//    latent statistics
scnt_int T[NUMMN]; 
stcnt_int *t[NUMMN];
stcnt_int t_data[DIM*NUMMN];
//  maxima used for building tables
int MAXN, MAXT;   

/*
 *  runtime statistics for reporting
 */
//  averages for a single run
int tcnt;
float Tave[NUMMN], *tave[NUMMN];
float tave_data[DIM*NUMMN];
//   runtime stats for bpar
float bave;
int   bcnt;
//   runtime stats for apar
float aave;
int   acnt;


/*==================================================
 * utilities
 *==================================================*/

static int sampleH() {
  double val = rng_unit(rng);
  int i;
  for (i=0; i<DIM; i++)
    if ( (val-=H[i])<=0 )
      break;
  if ( i>= DIM )
    i = DIM-1;
  return i;
}


/*==================================================
 *  sampling main
 *==================================================*/

void usage(int burnin, int ITER, int useN) {
  fprintf(stderr,"Commandline:  OPTION+\n");
  fprintf(stderr,
          "  OPTION is choice of:\n"
	  "   -a val,init  # discount par (default=%f)\n"
	  "   -b val,init  # concentration par (default=%f)\n"
	  "   -B burnin    # burnin before recording (default=25%%)\n"
	  "                #    end with 'ms' to set using millisec\n"
	  "   -C cycles    # recording cycles (default=%d)\n"
	  "   -H cycles    # cycles for sampling b (default=none)\n"
	  "   -I cycles    # cycles for sampling a (default=none)\n"
	  "   -N samples   # number of samples (default=%d)\n"
	  "   -s seed      # for random sampler\n"
	  "   -T max       # maximum t for building tables\n"
	  "   -v           # verbose\n",
	  apar, bpar, ITER, useN);
}

int main(int argc, char* argv[])
{
  int i, j, c, iter, ITER=200;
  unsigned long int seed=0;
  int bcycle = 0;
  float bstart = 0;
  int acycle = 0;
  float astart = 0;
  int burnin = 0;
  stable_t *ST = NULL;
  int useN = DIM*2;

  MAXN = 1;
  MAXT = MAXSTAB;
  /*
   *  default values for args
   */
  while ( (c=getopt(argc, argv,"a:b:B:C:I:hH:I:N:P:S:s:T:v"))>=0 ) {
    switch ( c ) {
    case 'h':
      usage(burnin?burnin:ITER/2, ITER, useN);
      exit(0);
    case 'b':
      if ( !optarg || sscanf(optarg,"%f,%f",&bpar, &bstart)<1 )
	yaps_quit("Need a valid 'b' argument\n");
      break;
    case 'a':
      if ( !optarg || sscanf(optarg,"%f,%f",&apar,&astart)<1 )
	yaps_quit("Need a valid 'a' argument\n");
      break;
    case 'H':
      if ( !optarg || sscanf(optarg,"%d",&bcycle)!=1 )
	yaps_quit("Need a valid 'G' argument\n");
      break;
    case 'T':
      if ( !optarg || sscanf(optarg,"%d",&MAXT)!=1 )
	yaps_quit("Need a valid 'T' argument\n");
      break;
     case 'I':
      if ( !optarg || sscanf(optarg,"%d",&acycle)!=1 )
	yaps_quit("Need a valid 'H' argument\n");
      break;
    case 'N':
      if ( !optarg || sscanf(optarg,"%d",&useN)!=1 )
	yaps_quit("Need a valid 'N' argument\n");
      break;
    case 'C':
      if ( !optarg || sscanf(optarg,"%d",&ITER)!=1 )
	yaps_quit("Need a valid 'C' argument\n");
       break;
    case 'B':
      if ( !optarg || sscanf(optarg,"%d",&burnin)!=1 )
	yaps_quit("Need a valid 'B' argument\n");
      break;
    case 's':
      if ( !optarg || sscanf(optarg,"%lu",&seed)!=1 )
	yaps_quit("Need a valid 's' argument\n");
      break;
    case 'v':
      verbose++;
      break;
#ifdef S_USE_THREADS
      case 'P':
      if ( !optarg || sscanf(optarg,"%u",&threads)!=1 )
	yaps_quit("Need a valid 'P' argument\n");
      break;
#endif
    default:
      yaps_message("Bad command line argument\n\n");
      usage(burnin?burnin:ITER/2, ITER, useN);
      exit(0);
    }
  }
  if ( useN>=MAXDATA ) 
    yaps_quit("N too large\n");
  
  if ( burnin==0 )
    burnin = ITER/2;
  else
    if ( burnin>=ITER-1 )
      yaps_quit("Burnin %d too large for cycles %d\n", burnin, ITER);

  yaps_message("Configuration details\n");
  yaps_message("=====================\n");
  /*
   *   set random number generator
   */
  if ( seed ) {
    rng_seed(rng,seed);
  } else {
    rng_time(rng,&seed);
  }
  yaps_message("Setting seed for data = %lu\n", seed);
    
  if ( acycle && apar==0 )
    apar = 0.5;
 
  yaps_message("Setting a=%f, b=%f, N=%d, D=%d\n", 
	       apar, bpar, useN, NUMMN);
  yaps_message("             burnin=%d,", burnin);
  yaps_message(" cycles=%d\n", ITER);
  
  /*
   *    fix pointers
   */
  for (j=0; j<NUMMN; j++) {
    n[j] = &n_data[j*DIM];
    t[j] = &t_data[j*DIM];
    tave[j] = &tave_data[j*DIM];
  }

  /*
   *    initialise everything
   */
  for (j=0; j<NUMMN; j++) {
    N[j] = useN;
    T[j] = 0;
    Tave[j] = 0;
    for (i=0; i<DIM; i++) {
      n[j][i] = 0;
      t[j][i] = 0;
      tave[j][i] = 0;
    }
  }

  /*
   *  fix base distribution, uniform
   */ 
  {
    for (i=0; i<DIM; i++) {
      H[i] = 1.0/DIM;
    }
  } 
  
  /*
   *     create data using a CRP to get initialisation for n[]
   */
  c = 0;
  for (j=0; j<NUMMN; j++) {
    int cc;
    i = sampleH();
    data[c++] = i;
    //  first entry always adds a table
    n[j][i]++;  
    t[j][i]++;
    T[j]++;
    for (cc=1; cc<N[j]; cc++) {
      float val = (cc+bpar)*rng_unit(rng);
      val -=  T[j]*apar+bpar;
      if ( val<=0 ) {
	//  new table
	i = sampleH();
	t[j][i]++;
	T[j]++;
      } else {
	for (i=0; i<DIM; i++) {
	  val -= n[j][i] - t[j][i]*apar;
	  if ( val<0 )
	    break;
	}
      }
      assert(i<DIM);
      n[j][i]++;  
      data[c++] = i;
    }
  }
  
  binit = bpar;
  
  /*
   *   record maximum entries in data
   *   do this where possible so that one can get the table
   *   sizes right
   *
   */
  MAXN = n[0][0]+1;
  MAXT = 1;
  for (j=0; j<NUMMN; j++) {
    for (i=0; i<DIM; i++) {
      if ( MAXN<=n[j][i] ) MAXN = n[j][i]+1;
      if ( MAXT<t[j][i] ) MAXT = t[j][i]*1.1+1;
    }
  }
  if ( MAXT>MAXN )
    MAXT = MAXN;
  
  yaps_message("Making S for N=%d M=%d a=%lf\n", MAXN,MAXT,apar);
  ST = S_make(MAXN, MAXT, MAXN, MAXTAB,
	      apar, S_STABLE | S_UVTABLE);
  if ( ST==NULL )
    yaps_quit("Making S failed!\n");
  S_report(ST,stdout);

  /*
   *    the seed only sets the data/sample,
   *    the seed for the simulation/Gibbs is always random
   */
  rng_free(rng);
  rng_time(rng,&seed);
  //yaps_message("Resetting seed = %lu\n", seed);
  
  /*
   *   report on initial data statistics
   */
  yaps_message("\nData sampled\n");
  yaps_message("============\n");
  for (j=0; j<NUMMN; j++) {
    yaps_message("n[%d] =", j);
    for (i=0; i<DIM; i++)
      yaps_message(" %d", n[j][i]);
    yaps_message(" = %d\n", N[j]);
    yaps_message("t[%d] =",j);
    for (i=0; i<DIM; i++)
      yaps_message(" %d", t[j][i]);
    yaps_message(" = %d\n", T[j]);
  }

  /*
   *   set the hyperparameters used in Gibbs,
   *   can be different to data
   */
  if ( bstart==0 ) 
    bstart = bpar;
  if ( astart==0 ) 
    astart = apar;

  //  initialise latent stats and reporting info
  for (j=0; j<NUMMN; j++) {
    T[j] = 0;
    Tave[j] = 0;
  }
  tcnt = 0;
  bave = 0;
  bcnt = 0;
  aave = 0;
  acnt = 0;

  bpar = bstart;
  if ( verbose && bcycle!=0 )
    yaps_message("Starting with initial b=%f\n", bpar);
  
  apar = astart;
  if ( verbose && acycle!=0 )
    yaps_message("Starting with initial a=%f\n", apar);
  
  for (j=0; j<NUMMN; j++) {
    for (i=0; i<DIM; i++) {
      tave[j][i] = 0;
      t[j][i] = 0;
      if ( n[j][i]>0 ) {
	/*
	 *  initialise to a single table
	 */
	t[j][i] = 1;
	T[j]++;
      }
    }
  }

  for ( iter=0; iter<ITER; iter++) {
    /*
     *   sampling with table indicators
     */
    c = 0;
    for (j=0; j<NUMMN; j++) {
      int cc;
      for (cc=0; cc<N[j]; cc++) {
	float one;
	i = data[c++];
	assert(n[j][i]);
	if ( n[j][i]==1 )
	  //    this indicator must always be 1, no sampling
	  continue;
	//   sample whether it contributes to a table
	if ( t[j][i]>1 &&
	     (n[j][i]-1)*rng_unit(rng)<(t[j][i]-1) ) {
	  t[j][i]--;  
	  T[j]--;
	} 
	assert(t[j][i]<n[j][i]);
	//  sample new table indicator
	one = H[i] * (bpar + T[j]*apar) * (t[j][i]) / (n[j][i]-t[j][i]+1)
	  * S_V(ST, n[j][i],t[j][i]+1);
	if ( rng_unit(rng) < one/(one+1.0) ) {
	  t[j][i]++;
	  T[j]++;
	} 
      }
    }
    
    /*
     *   one major cycle of Gibbs sampler finished
     */
    if ( verbose>1 ) {
      for (j=0; j<NUMMN; j++) {
	for (i=0; i<DIM; i++)
	  printf(" %d", t[j][i]);
	printf(" = %d\n", T[j]);
      }
    }
    /*
     *     sample & record b 
     */
    if ( bcycle!=0 && iter%bcycle==0 ) {
      //   Gibbs on bpar (concentration par) too
      if ( bcycle<0 ) {
	  int bc = -bcycle;
	  for (bc-- ; bc>0; bc--)
	    bpar = sampleb(bpar, 1, PB_shape, PB_scale, N, T, 
			   apar, rng, 1, 1);
      }
      bpar = sampleb(bpar, 1, PB_shape, PB_scale, N, T, 
		     apar, rng, 1, 1);
      if ( iter>=burnin ) {
	bave += bpar;
	bcnt ++;
      }
    }
    /*
     *     sample & record a
     */
    if ( acycle!=0 && iter%acycle==0 ) {
      int dimI[NUMMN];
      double dimb[NUMMN];
      for (j=0; j<NUMMN; j++)  {
	dimI[j] = DIM;
	dimb[j] = bpar;
      }
      //   Gibbs on apar (discount par) too
      if ( acycle<0 ) {
	int bc = -acycle;
	for (bc-- ; bc>0; bc--)
	  apar = samplea(apar, NUMMN, dimI, T, n, t, NULL, dimb, rng, 1, 1);
      }
      apar = samplea(apar, NUMMN, dimI, T, n, t, NULL, dimb, rng, 1, 1);
      if ( iter>=burnin ) {
	aave += apar;
	acnt ++;
      }
      if ( verbose>1 )
	yaps_message("Extending S for a=%lf\n", apar);
      if ( S_remake(ST,apar) )
	yaps_message("Extending S failed\n");
    }
    /*
     *    full statistics collection
     */
    if ( iter>=burnin ) {
      for (j=0; j<NUMMN; j++) {
	for (i=0; i<DIM; i++) {
	  tave[j][i] += t[j][i];
	}
	Tave[j] += T[j];
      }
      tcnt ++;
    }
  }	  
  
  /*
   *     report for this experiment
   */
  yaps_message("\nEstimates\n");
  yaps_message("=========\n");
  for (j=0; j<NUMMN; j++) {
    yaps_message("t[%d] = ", j);
    for (i=0; i<DIM; i++)
      yaps_message(" %.2f", tave[j][i]/tcnt);
    yaps_message("\nT[%d]=%.2f\n", j, Tave[j]/tcnt);
  }
  if ( bcycle!=0 && bcnt>0 ) 
    yaps_message("\nb=%.2f", bave/bcnt);
  if ( acycle!=0 && acnt>0 ) 
    yaps_message("\na=%.3f", aave/acnt);
  yaps_message("\n");
  
  S_free(ST);
  rng_free(rng);
  return 0;
}
