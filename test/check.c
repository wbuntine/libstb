/*
 *  WARNING:  ignore this file, use "demo.c" instead
 *
 * Extensive Dirichlet/Pitman-Yor Process sampling test
 * Copyright (C) 2011-2012 Wray Buntine and Lan Du
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
 * This is an extensive suite (demo.c is a cut down version of this)
 * for evaluating estimation of multinomial/Pitman-Yor Process sampling.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif


#include "yaps.h"
#include "stable.h"
#include "lgamma.h"
#include "digamma.h"
#include "psample.h"
/*
 *    the RNG is GSL (usually) and all done through rng_*() interface
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

int SliceSimple(double *xp, double (*post)(double, void *), 
		double *bounds,	rngp_t rng, int loops, void *pars);


/*==================================================
 * global variables
 *
 *     This is a very primitive setup to test a single
 *     Pitman-Yor.
 *==================================================*/

rngp_t rng;
int verbose = 0;

/*
 *    by default, do slice sampling for b,
 *    otherwise use ARS
 *       NB.  they seem to be similar in performance after basic tuning
 */
int use_ars = 0;
/*
 *    dimensions and constants
 */
//    H() is discrete distribution with this dimension
#define DIM 5
//    array limits
#define MAXSTAB 10000
#define MAXTAB 10000
#define MAXDATA 100000

/*
 *   base distribution
 */
enum BaseType { BaseUniform, BaseLinear, BaseSlowLinear, BaseDirichlet } 
  Htype = BaseUniform;

/*
 *    SampleSA = sampling for seating arrangements
 *    SampleHSA = histogrammed sampling for seating arrangements
 *    SampleCT = collapsed table sampler
 *    SampleCTW = collapsed table sampler wih windowing
 *    SampleTI = table indicator sampler
 */
enum SampleType { SampleSA, SampleCT, SampleCTW, SampleTI, SampleHSA };

/*
 *   parameters and hyper parameters
 */
float  apar = 0.0, bpar = 1.0, binit;
#define PB_shape 1.1
#define PB_scale DIM
//    base distribution
float H[DIM];
//    window size for bounded CT sampler
int TWINDOW=10;

/*
 *   data
 */
int data[MAXDATA];
scnt_int N=DIM*5;
//    statistics
scnt_int n[DIM];
//    latent statistics
scnt_int T=DIM;
stcnt_int t[DIM], *(m[DIM]), one[DIM];
stcnt_int mm[DIM*MAXTAB];   // storage space for m[]
//    first storage
int f[DIM];
int MAXN, MAXT;       //  maxima used for building tables

/*
 *  runtime statistics for reporting
 */
//  averages for a single run
int tcnt;
float Tave, tave[DIM];
//  averages over runs
int   avecnt;
float Taveave, taveave[DIM];

#define REPAVE 3
float *(runave[REPAVE]);
float *(sqrunave[REPAVE]);

//   runtime stats for bpar
float bave, baveave;
int   bcnt;
//   runtime stats for apar
float aave, aaveave;
int   acnt;


/*==================================================
 * utilities
 *==================================================*/

static int my_gettime(struct timespec *ts) {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
   clock_serv_t cclock;
   mach_timespec_t mts;
   host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
   clock_get_time(cclock, &mts);
   mach_port_deallocate(mach_task_self(), cclock);
   ts->tv_sec = mts.tv_sec;
   ts->tv_nsec = mts.tv_nsec;
   return 0;
#else
   return clock_gettime(CLOCK_REALTIME, ts);
#endif
}

static int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p)
{
  return ((timeA_p->tv_sec-timeB_p->tv_sec) * 1000000000) 
    + (timeA_p->tv_nsec - timeB_p->tv_nsec);
}

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
 *  sampling a
 *  prior on a is uniform
 *
 *  either sample on the simple Gamma-based posterior,
 *  or switch to the library function
 *==================================================*/

#define LGCACHE
static double aterms(double x, void *mydata) {
  int i, tt;
  double val = T * log(x) + lgamma(T+bpar/x) - lgamma(bpar/x);
#ifdef LGCACHE
  struct gcache_s lgp;
  gcache_init(&lgp, x);
#else
  double lg1x = lgamma(1-x);
#endif
  for (i=0; i<DIM; i++)
    if ( n[i]>1 )
      for (tt=0; tt<t[i]; tt++) 
	if ( m[i][tt]>1 ) {
#ifdef LGCACHE
	  val += gcache_value(&lgp, m[i][tt]);
#else
	  val += lgamma(m[i][tt]-x) - lg1x;
#endif
	}
  return val;
}

/*
 *    use  prior a is uniform
 */
#define SQUEEZEA 0.2
float mysamplea(float mya, enum SampleType sampler, stable_t *ST ) {
  double result;
  double inita[3] = {A_MIN,1,A_MAX};
  int i;
  if ( (sampler==SampleSA||sampler==SampleHSA) ) {
    /*
     *   these samplers use the simple posterior,
     *   no Stirling numbers
     */
    inita[1] = mya;
    if ( fabs(inita[1]-A_MAX)/A_MAX<0.00001 ) {
      inita[1] = A_MAX*0.999 + A_MIN*0.001;
    }
    if ( fabs(inita[1]-A_MIN)/A_MIN<0.00001 ) {
      inita[1] = A_MIN*0.999 + A_MAX*0.001;
    }
#ifdef SQUEEZEA
    /*
     *  bound current move to be less than SQEEZEA
     */
    if ( inita[1]-SQUEEZEA>A_MIN )  inita[0] = inita[1]-SQUEEZEA;
    if ( inita[1]+SQUEEZEA<A_MAX )  inita[2] = inita[1]+SQUEEZEA;
#endif
    //  yaps_message("Sample a ~ ARS for a=[%lf,%lf,%lf]\n", 
    //	      inita[0],inita[1],inita[2]);
    MAXT = t[0];
    for (i=1; i<DIM; i++)
      if ( t[i]>MAXT )
	MAXT = t[i]+1;
    if ( use_ars ) {
      arms_simple(3,  inita, inita+2, aterms, NULL, 0, inita+1, &result);
      if ( result<inita[0] ||  result>inita[2] )
	yaps_quit("Arms_simple(apar) returned value out of bounds\n");
    } else {
      result = mya;
      inita[1] = A_MAX;
      if ( SliceSimple(&result, aterms, inita, rng, 1, NULL) )
	yaps_quit("SliceSimple error\n");
    }
    mya = result;
    if ( verbose>1 )
      yaps_message("Sample a ~ A(%d) = %f\n", T, mya);
    return mya;
  } else {
    scnt_int *np = &n[0];
    stcnt_int *nt = &t[0];
    double myb = bpar;
    int myK = DIM;
#ifdef SAMPLEA_M
    return samplea2(mya, ST, 1, &myK, &T, &np, &nt, NULL, &myb, rng, 1, 1);
#else
    return samplea(mya, 1, &myK, &T, &np, &nt, NULL, &myb, rng, 1, 1);
#endif
  }
}


/*==================================================
 *  sampling main
 *==================================================*/

void usage(int burnin, int ITER) {
  fprintf(stderr,"Commandline:  OPTION+\n");
  fprintf(stderr,
          "  OPTION is choice of:\n"
	  "   -A           # use ARS, not slice sampling for a, b\n"
	  "   -a val,init  # discount par (default=%f)\n"
	  "   -b val,init  # concentration par (default=%f)\n"
	  "   -B burnin    # burnin before recording (default=25%%)\n"
	  "                #    end with 'ms' to set using millisec\n"
	  "   -c maxrel    # bound t sampling when ratio drops\n"
	  "   -C cycles    # recording cycles (default=%d)\n"
	  "                #    end with 'ms' to set using millisec\n"
	  "   -H cycles    # cycles for sampling b (default=none)\n"
	  "   -I cycles    # cycles for sampling a (default=none)\n"
	  "   -N samples   # number of samples (default=%d)\n"
	  "   -p repcyc    # print current mean t per repcyc cycles\n"
	  "   -r repeat    # times to repeat wth different sample\n"
	  "   -s seed      # for random sampler\n"
	  "   -S Type      # Type=CT (collapsed table sampler)\n"
	  "                #     =CTW (windowed collapsed table sampler)\n"
	  "                #     =TI (table indicator sampler)\n"
	  "                #     =SA (seating arrangement sampler)\n"
	  "                #     =HSA (histogrammed SA sampler)\n"
	  "   -T max       # maximum t for building tables\n"
	  "   -v           # verbose\n"
	  "   -w size      # window of [-size,size] sampled for -SCTW (def=%d)\n",
	  apar, bpar, ITER, N, TWINDOW);
}

int main(int argc, char* argv[])
{
  int i, c, iter, ITER=200;
  unsigned long seed=0;
  int ITERTIME = 0;
  int printmean = 0;
  int bcycle = 0;
  float bstart = 0;
  int acycle = 0;
  float maxrel = 1e30;
  float astart = 0;
  int REDO = 1;
  int redo = 1;
  enum SampleType sampler = SampleCT; 
  int burnin = 0;
  int burnintime = 0;
  stable_t *ST = NULL;
  struct timespec tstart, tend;
  /*
   *   4 values kept for the first run only
   *      time, mean-T, mean-b, mean-1
   */
#define REPVALS 4
  double *repval=NULL;

  MAXN = 1;
  MAXT = MAXSTAB;
  /*
   *  default values for args
   */
  while ( (c=getopt(argc, argv,"Aa:b:B:c:C:I:hH:I:p:N:q:r:S:s:T:vw:"))>=0 ) {
    switch ( c ) {
    case 'c':
     if ( !optarg || sscanf(optarg,"%f",&maxrel)!=1 )
	yaps_quit("Need a valid 'c' argument\n");
      break;
    case 'p':
     if ( !optarg || sscanf(optarg,"%d",&printmean)!=1 )
	yaps_quit("Need a valid 'p' argument\n");
      break;
    case 'h':
      usage(burnin?burnin:ITER/2, ITER);
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
	yaps_quit("Need a valid 'H' argument\n");
      break;
    case 'T':
      if ( !optarg || sscanf(optarg,"%d",&MAXT)!=1 )
	yaps_quit("Need a valid 'T' argument\n");
      break;
     case 'I':
      if ( !optarg || sscanf(optarg,"%d",&acycle)!=1 )
	yaps_quit("Need a valid 'I' argument\n");
      break;
    case 'r':
      if ( !optarg || sscanf(optarg,"%d",&REDO)!=1 )
	yaps_quit("Need a valid 'r' argument\n");
      break;
    case 'N':
      if ( !optarg || sscanf(optarg,"%d",&N)!=1 )
	yaps_quit("Need a valid 'N' argument\n");
      break;
    case 'C':
      if ( !optarg || sscanf(optarg,"%d",&ITER)!=1 )
	yaps_quit("Need a valid 'C' argument\n");
      if ( strcmp(optarg+strlen(optarg)-2,"ms")==0 ) {
	ITERTIME = ITER;
	ITER = 100000000;
      }
       break;
    case 'B':
      if ( !optarg || sscanf(optarg,"%d",&burnin)!=1 )
	yaps_quit("Need a valid 'B' argument\n");
      if ( strcmp(optarg+strlen(optarg)-2,"ms")==0 ) {
	burnintime = burnin;
	burnin = 0;
      }
      break;
    case 's':
      if ( !optarg || sscanf(optarg,"%lu",&seed)!=1 )
	yaps_quit("Need a valid 's' argument\n");
      break;
    case 'w':
      if ( !optarg || sscanf(optarg,"%d",&TWINDOW)!=1 )
	yaps_quit("Need a valid 'w' argument\n");
      break;
    case 'S':
      if ( !optarg ) 
	yaps_quit("Need a valid 'S' argument\n");
      if ( strcmp(optarg,"CTW")==0 )
	sampler = SampleCTW;
      else if ( strcmp(optarg,"CT")==0 )
	sampler = SampleCT;
      else if ( strcmp(optarg,"TI")==0 )
	sampler = SampleTI;
      else if ( strcmp(optarg,"SA")==0 )
	sampler = SampleSA;
      else if ( strcmp(optarg,"HSA")==0 )
	sampler = SampleHSA;
      else
	yaps_quit("Need a valid 'S' argument\n");
      break;
    case 'A':
      use_ars++;
      break;
    case 'v':
      verbose++;
      break;
    default:
      yaps_message("Bad command line argument\n\n");
      usage(burnin?burnin:ITER/2, ITER);
      exit(0);
    }
  }
  if ( N>=MAXDATA ) 
    yaps_quit("N too large\n");
  
  /*
   *   print an initial report about the
   *   configuration used
   */
  switch ( sampler ) {
  case SampleSA:
    yaps_message("Seating arrangements sampler\n");
    break;
  case SampleHSA:
    yaps_message("Tabulated seating arrangements sampler\n");
    break;
  case SampleTI:
    yaps_message("Table indicator sampler\n");
    break;
  case SampleCT:
    yaps_message("Collapsed table sampler\n");		
    break;
  case SampleCTW:
    yaps_message("Collapsed table sampler with windowing\n");		
    break;
  default:
    yaps_quit("Bad sampler\n");
  }
  switch ( Htype ) {
  case BaseUniform:
    yaps_message("Uniform base distribution of dimension K=%d\n", DIM);		
    break;
  case BaseLinear:
    yaps_message("Linear (ramp down) base distribution of dimension K=%d\n", DIM);		
    break;
  case BaseSlowLinear:
    yaps_message("Slow linear (ramp down) base distribution of dimension K=%d\n", DIM);		
    break;
  case BaseDirichlet:
    yaps_message("Truncated stick breaking distribution of dimension K=%d\n", DIM);		
    break;
    
  }
  
  /*
   *  horrible fiddle ... adjust iterations
   *  to make the different reporting/sampling "even"
   */
  if ( bcycle!=0 ) 
    while ( (ITER%(bcycle))!=0 ) ITER++;

  if ( ITERTIME>0 && burnin==0 ) {
    burnintime = ITERTIME/4;
  }
  if ( burnintime==0 ) {
    if ( burnin==0 )
      burnin = ITER/2;
    else
      if ( burnin>=ITER-1 )
	yaps_quit("Burnin %d too large for cycles %d\n", burnin, ITER);
  }

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
 
  yaps_message("Running with a=%f, b=%f, N=%d\n", apar, bpar, N);
  if ( burnintime ) 
    yaps_message("             burnin=%dms,",
		burnintime);
  else
    yaps_message("             burnin=%d,", 
		burnin);
  if ( ITERTIME ) 
    yaps_message(" cycles=%dms", ITERTIME);
  else
    yaps_message(" cycles=%d", ITER);
  yaps_message(", repeats=%d\n", REDO);
  
  /*
   *    initialise everything
   */
  T = 0;
  Taveave = 0;
  baveave = 0;
  aaveave = 0;
  Tave = 0;
  avecnt = REDO;
  if ( printmean ) {
    int REP=(ITER-burnin)/printmean;
    int ii;
    for (ii=0; ii<REPAVE; ii++) {
      runave[ii] = calloc(REP,sizeof(*runave[0]));
      sqrunave[ii] = calloc(REP,sizeof(*sqrunave[0]));
    }
    repval = calloc(REPVALS*REP,sizeof(*repval));
  }
  for (i=0; i<DIM; i++) {
    n[i] = 0;
    t[i] = 0;
    f[i] = -1;
    m[i] = &mm[i*MAXTAB];
    m[i][0] = 0;
    tave[i] = 0;
    taveave[i] = 0;
  }

  /*
   *  fix base distribution
   */ 
  if ( Htype!=BaseDirichlet ) {
    /*
     *  uniform or linear ramping down
     */
    double tot = 0;
    for (i=0; i<DIM; i++) {
      if ( Htype==BaseLinear )
	tot += H[i] = DIM-i;
      else if ( Htype==BaseSlowLinear )
	tot += H[i] = DIM*1.5-i;
      else
	tot += H[i] = 1.0;
    }
    for (i=0; i<DIM; i++) 
      H[i] /= tot;
  } else {
    /*
     *  truncated stick breaking distribution
     */
    double left = 1;
    for (i=0; i<DIM; i++) {
      double val = rng_beta(rng, 1, bpar);
      H[i] = left*val;
      left *= (1.0-val);
    }
    H[DIM-1] += left;
  }

  /*
   *     create data using a CRP to get initialisation for n[]
   */
  i = sampleH();
  data[0] = i;
  f[i] = 0;
  //  first entry always adds a table
  n[i]++;  
  t[i]++;
  T++;
  for (c=1; c<N; c++) {
    float val = (c+bpar)*rng_unit(rng);
    val -=  T*apar+bpar;
    if ( val<=0 ) {
      //  new table
      i = sampleH();
      if ( f[i]<0 )
	f[i] = c;
      t[i]++;
      T++;
    } else {
      for (i=0; i<DIM; i++) {
	val -= n[i] - t[i]*apar;
	if ( val<0 )
	  break;
      }
    }
    assert(i<DIM);
    n[i]++;  
    data[c] = i;
  }
  
  binit = bpar;
  
  /*
   *   record maximum entries in data
   */
  MAXN = n[0]+1;
  for (i=1; i<DIM; i++) {
    if ( MAXN<=n[i] ) MAXN = n[i]+1;
    if ( MAXT<t[i] ) MAXT = t[i]*1.1+1;
  }
  if ( MAXT>MAXN )
    MAXT = MAXN;
  
  yaps_message("Making S for N=%d M=%d a=%lf\n", MAXN,MAXT,apar);
  ST = S_make(MAXN,MAXT,MAXN,MAXTAB,
	      apar, S_STABLE | (sampler==SampleTI?S_UVTABLE:0));
  if ( ST==NULL )
    yaps_quit("Making S failed!\n");
  S_report(ST,stderr);

  /*
   *    the seed only sets the data/sample,
   *    the seed for the simulation/Gibbs is always random
   */
  rng_free(rng);
  rng_time(rng,&seed);
  yaps_message("Resetting seed = %lu\n", seed);
  
  /*
   *   report on initial data statistics
   */
  yaps_message("n[] =");
  for (i=0; i<DIM; i++)
    yaps_message(" %d", n[i]);
  yaps_message(" = %d\n", N);
  yaps_message("t[] =");
  for (i=0; i<DIM; i++)
    yaps_message(" %d", t[i]);
  yaps_message(" = %d\n", T);

  /*
   *   set the hyperparameters used in Gibbs,
   *   can be different to data
   */
  if ( bstart==0 ) 
    bstart = bpar;
  if ( astart==0 ) 
    astart = apar;

  /*
   *     we repeat the experiment REDO times
   */
  if ( my_gettime(&tstart)!=0 )
    yaps_sysquit("Timer failed\n");
  for (redo=REDO ; redo>0; redo--) {
    int repi = 0;
    //  initialise latent stats and reporting info
    T = 0;
    Tave = 0;
    tcnt = 0;
    bave = 0;
    bcnt = 0;
    aave = 0;
    acnt = 0;
    /*
     *   initialise bpar to a range of values if sampling b
     */
    if ( bcycle!=0 && REDO>1 ) 
      bpar = bstart/3.0  +  bstart*(3.0 - 1/3.0) * (REDO-redo)/REDO;
    else
      bpar = bstart;
    if ( verbose && bcycle!=0 )
      yaps_message("Restarting with initial b=%f\n", bpar);

    /*
     *   initialise apar to a range of values if sampling a
     */
    if ( acycle!=0 && REDO>1 ) 
      apar = A_MIN + (A_MAX-A_MIN)*rng_unit(rng);
    else
      apar = astart;
    if ( verbose && acycle!=0 )
      yaps_message("Restarting with initial a=%f\n", apar);

    for (i=0; i<DIM; i++) {
      t[i] = 0;
      if ( n[i]>0 ) {
	/*
	 *  initialise to a single table
	 */
	t[i] = 1;
	tave[i] = 0;
	T++;
      }
      one[i] = 0;
      m[i][0] = 0;
      if ( n[i]>0 ) {
	if ( sampler==SampleSA || (sampler==SampleHSA && n[i]>1 )) {
	  m[i][0] = n[i];
	  m[i][1] = 0;
	} else if ( sampler==SampleHSA && n[i]==1 ) {
	  one[i] = 1;
	} 
      }
    }
      
    for ( iter=0; iter<ITER; iter++) {
      if ( sampler == SampleSA ) {
	/*
	 *   standard sampling seating arrangements
	 *   where table sizes are explicitly kept
	 */
	for (c=0; c<N; c++) {
	  int tt, nn;
	  float val;
	  i = data[c];
	  //   sample table
	  assert(n[i]);
	  nn = n[i]*rng_unit(rng);
	  for (tt=0; tt<t[i]; tt++) {
	    if ( (nn-=m[i][tt])<0 )
	      break;
	  }
	  if ( verbose>1 )
	    yaps_message("i=%d withdraw t=%d/%d\n", i, tt, t[i]);
	  assert(tt<t[i]);
	  //   remove from table
	  m[i][tt]--;
	  assert( m[i][tt]>=0);
	  if ( m[i][tt]==0 ) {
	    t[i]--;
	    m[i][tt] = m[i][t[i]];  // shuffle
	    T--;
	  }
	  //  sample new table, noting "i" is removed so (N-1) etc.
	  val = (T*apar+bpar)*H[i] + n[i]-1 - t[i]*apar;
	  if ( verbose>1 )
	    yaps_message("sampling  new=%lf  vs  old=%lf\n",
			(T*apar+bpar)*H[i],
			n[i]-1 - t[i]*apar);
	  val *= rng_unit(rng);
	  for (tt=0;  tt<t[i] && (val-=m[i][tt]-apar)>0; tt++) ;
	  if ( tt>=t[i] ) {
	    if ( tt>=MAXTAB )
	      yaps_quit("Too many tables during sampling\n");
	    m[i][t[i]] = 1;
	    t[i]++;
	    T++;
	    if ( verbose>1 )
	      yaps_message("i=%d new table t=%d\n", i, t[i]); 
	  } else {
	    assert(tt<t[i]);
	    m[i][tt]++;	
	    if ( verbose>1 )
	      yaps_message("i=%d increment t=%d/%d\n", i, tt, t[i]); 
	  }
	}
      } else if ( sampler == SampleHSA ) {
	/*
	 *   standard sampling seating arrangements
	 *   where table sizes are explicitly kept;
	 *   but the tables with size 1 are histogrammed
	 */
	for (c=0; c<N; c++) {
	  int tt;
	  double val;
	  i = data[c];
	  {
	    int nt = one[i];
	    for (tt=t[i]-one[i]-1;  tt>=0; tt--) 
	      nt += m[i][tt];
	    assert(nt==n[i]);
	  }
	  //   sample table
	  val = n[i]*rng_unit(rng);
	  val -= one[i];
	  if ( one[i]>0 && val<0 ) {
	    t[i]--;
	    one[i]--;
	    T--;
	  } else {
	    for (tt=t[i]-one[i]-1; tt>=0; tt--) {
	      if ( (val-=m[i][tt])<0 )
		break;
	    }
	    if ( tt<0 ) tt = 0;
	    //   remove from table
	    m[i][tt]--;
	    assert( m[i][tt]>=1);
	    if ( m[i][tt]==1 ) {
	      //  convert to a one
	      one[i]++;
	      m[i][tt] = m[i][t[i]-one[i]];  // shuffle out
	    }
	  }
	  {
	    int nt = one[i];
	    for (tt=t[i]-one[i]-1;  tt>=0; tt--) 
	      nt += m[i][tt];
	    assert(nt==n[i]-1);
	  }
	  //  sample new table, noting "i" is removed so (N-1) etc.
	  val = (T*apar+bpar)*H[i] + n[i]-1 - t[i]*apar;
	  val *= rng_unit(rng);
	  val -= one[i]*(1-apar);
	  if ( val<0 && one[i]>0 ) {
	    //  turn a one into a regular table
	    m[i][t[i]-one[i]] = 2;
	    one[i]--;
	  } else {
	    for (tt=t[i]-one[i]-1;  tt>=0 && (val-=m[i][tt]-apar)>0; tt--) ;
	    if ( tt<0 ) {
	      /*
	       *  new table, necessarily a one
	       */
	      t[i]++;
	      if ( t[i]>=MAXTAB )
		yaps_quit("Too many tables during sampling\n");
	      one[i]++;
	      T++;
	    } else {
	      //  increment a regular table
	      assert(tt<t[i]);
	      assert(m[i][tt]>1);
	      m[i][tt]++;	
	    }
	  }
	  {
	    int nt = one[i];
	    for (tt=t[i]-one[i]-1;  tt>=0; tt--) 
	      nt += m[i][tt];
	    assert(nt==n[i]);
	  }
	}
      } else if ( sampler == SampleTI ) {
	/*
	 *   sampling with table indicators
	 */
	for (c=0; c<N; c++) {
	  float one;
	  i = data[c];
	  assert(n[i]);
	  if ( n[i]==1 || c==f[i] )
	    //    this indicator must always be 1, no sampling
	    continue;
	  //   sample whether it contributes to a table
	  if ( t[i]>1 && (n[i]-1)*rng_unit(rng)<(t[i]-1) ) {
	    t[i]--;  
	    T--;
	  } 
	  assert(t[i]<n[i]);
	  //  sample new table indicator
	  one = H[i] * (bpar + T*apar) * (t[i]) / (n[i]-t[i]+1)
	    * S_V(ST, n[i],t[i]+1);
	  if ( rng_unit(rng) < one/(one+1.0) ) {
	    t[i]++;
	    T++;
	  } 
	}
      } else if ( sampler==SampleCT ) { 
	/*
	 *   collapsed sampler 
	 */
	for (i=0; i<DIM; i++) {
	  if ( n[i]>0 ) {
	    int tt;
	    double  try[MAXTAB];
	    double  hval = 1.0;
	    double  sval = S_S(ST,n[i],t[i]); 
	    double  tot;
	    int t_high=1;   // records highest value in loop
	    tot = try[1] = exp(S_S(ST,n[i],1)-sval);
	    for (tt=2; tt<=n[i] && tt<MAXTAB; tt++) {
	      hval *= H[i] * (bpar + (T-t[i]+tt-1)*apar);
	      try[tt] = hval * exp(S_S(ST,n[i],tt)-sval);
	      tot += try[tt];
	      if ( try[tt]> try[tt-1] ) 
		t_high=tt;
	      else 
		/*
		 *  early stopping of loop:
		 *        way after last t used
		 *        much smaller than largest value
		 */
		if ( tt>t[i]+3 && try[t_high]>maxrel*try[tt] )
		    break;
	    }
	    t_high = tt;  //  sets upper bound for sampling
	    tot *= rng_unit(rng);
	    tt = 1;
	    for (tt=1; tt<t_high && (tot-=try[tt])>0; tt++) ;
	    assert(tt<=n[i]);
	    T += tt-t[i];
	    t[i] = tt;
	  }
	}
      } else if ( sampler==SampleCTW ) { 
	/*
	 *   collapsed sampler with windowing
	 */
	for (i=0; i<DIM; i++) {
	  if ( n[i]>0 ) {
	    int tt;
	    double  try[MAXTAB];
	    double  hval;
	    double  sval = S_S(ST,n[i],t[i]); 
	    double  tot;
	    int maxt = MAXTAB-1, mint = 1;
	    if ( maxt>n[i] ) maxt = n[i];
	    if ( maxt>t[i]+TWINDOW ) maxt = t[i]+TWINDOW;
	    if ( mint<t[i]-TWINDOW ) mint = t[i]-TWINDOW;
	    tot = try[t[i]] = 1;
	    hval = 1.0;
	    for (tt=t[i]+1; tt<=maxt; tt++) {
	      hval *= H[i] * (bpar + (T-t[i]+tt-1)*apar);
	      try[tt] = hval * exp(S_S(ST,n[i],tt)-sval);
	      tot += try[tt];
	    }
	    hval = 1.0;
	    for (tt=t[i]-1; tt>=mint; tt--) {
	      hval /= H[i] * (bpar + (T-t[i]+tt)*apar);
	      try[tt] = hval * exp(S_S(ST,n[i],tt)-sval);
	      tot += try[tt];
	    }
	    tot *= rng_unit(rng);
	    for (tt=t[i]; tt<=maxt && (tot-=try[tt])>0; tt++) ;
	    if ( tot>0 )
	      for (tt=t[i]-1; tt>=mint && (tot-=try[tt])>0; tt--) ;
	    assert(tt<=n[i]);
	    assert(tt>=1);
	    T += tt-t[i];
	    t[i] = tt;
	  }
	}
      } 
      /*
       *   one major cycle of Gibbs sampler finished
       */
      if ( verbose>1 ) {
	for (i=0; i<DIM; i++)
	  printf(" %d", t[i]);
	printf(" = %d\n", T);
      }
      /*
       *     sample & record b 
       */
      if ( bcycle!=0 && iter%bcycle==0 ) {
	//   Gibbs on bpar (concentration par) too
	if ( bcycle<0 ) {
	  int bc = -bcycle;
	  for (bc-- ; bc>0; bc--)
	    bpar = sampleb(bpar, 1, PB_shape, PB_scale, &N, &T, 
			   apar, rng, 1, 1);
	}
	bpar = sampleb(bpar, 1, PB_shape, PB_scale, &N, &T, 
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
	//   Gibbs on apar (discount par) too
	if ( acycle<0 ) {
	  int bc = -acycle;
	  for (bc-- ; bc>0; bc--)
	    apar = mysamplea(apar,sampler,ST);
	}
	apar = mysamplea(apar,sampler,ST);
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
	for (i=0; i<DIM; i++) {
	  tave[i] += t[i];
	}
	tcnt ++;
	Tave += T;
	if ( printmean && (printmean==1 || iter%printmean==0) ) {
	  if ( redo==1 ) {
	    if ( tcnt==1 )
	      my_gettime(&tstart);
	    my_gettime(&tend);
	    repval[repi*REPVALS] = timespecDiff(&tend, &tstart)/1.0e6;
	    repval[repi*REPVALS+1] = Tave/tcnt;
	    if ( bcycle!=0 ) 
	      repval[repi*REPVALS+2] = bave/bcnt;
	    else
	      repval[repi*REPVALS+2] = bpar;
	    if ( acycle!=0 ) 
	      repval[repi*REPVALS+3] = aave/acnt;
	    else
	      repval[repi*REPVALS+3] = apar;
	  }
	  runave[0][repi] += Tave/tcnt;
	  sqrunave[0][repi] += (Tave/tcnt)*(Tave/tcnt);
	  if ( bcycle!=0 ) {
	    runave[1][repi] += bave/bcnt;
	    sqrunave[1][repi] += (bave/bcnt)*(bave/bcnt);
	  }
	  if ( acycle!=0 ) {
	    runave[2][repi] += aave/acnt;
	    sqrunave[2][repi] += (aave/acnt)*(aave/acnt);
	  }
	  repi++;
	}
      }
      if ( (ITERTIME>0 || burnintime>0) && redo==REDO && iter%10==0 ) {
	/*
	 * on the first run, we have to check stopping time
	 */
	double dt;
	my_gettime(&tend);
	dt = ((tend.tv_sec - tstart.tv_sec) * 1000.0
	      + (tend.tv_nsec - tstart.tv_nsec)/1000000.0);	
	if ( burnintime>0 && dt>burnintime ) {
	  yaps_message("Resetting burnin to %d for time bound\n",
		     iter+1);
	  /*
	   *   forces a stop now, a well a affecting later runs
	   */
	  burnin = iter+1;
	  burnintime = 0;
	}	
	if ( ITERTIME>0 && dt>ITERTIME ) {
	  yaps_message("Resetting iterations to %d for time bound\n",
		     iter+1);
	  /*
	   *   forces a stop now, a well a affecting later runs
	   */
	  ITER = iter+1;
	  ITERTIME = 0;
	}
      }
    }	  
    /*
     *   end of one sampling run;
     *   now need to record details for multiple sampling runs
     */

    //  record runtime stats
    for (i=0; i<DIM; i++) {
      tave[i] /= tcnt;
      taveave[i] += tave[i];
    }
    Tave /= tcnt;
    Taveave += Tave;
    if ( bcycle!=0 && bcnt>0 ) 
      baveave += bave/bcnt;
    if ( acycle!=0 && acnt>0 ) 
      aaveave += aave/acnt;

    if ( verbose ) {
      //   report for this experiment
      yaps_message("Ave: ");
      if ( (DIM>20 && verbose) || REDO==1 ) 
	for (i=0; i<DIM; i++)
	  yaps_message(" %f", tave[i]);
      yaps_message(" T=%f", Tave);
      if ( bcycle!=0 && bcnt>0 ) 
	yaps_message(", b=%f", bave/bcnt);
      if ( acycle!=0 && acnt>0 ) 
	yaps_message(", a=%f", aave/acnt);
      yaps_message("\n");
    }
  }
  /*
   *   end of sampling, now do reports
   */
  if ( my_gettime(&tend)!=0 )
    yaps_sysquit("Timer failed\n");

  {
    uint64_t timeElapsed = timespecDiff(&tend, &tstart);
    yaps_message("\nTime = %.1lf ms\n", timeElapsed/1.0e6);

  }
  if ( printmean ) {
    int REP=(ITER-burnin)/printmean;
    for (i=0; i<REP; i++) {
      float stderr;
      float mymean;
      runave[0][i] /= REDO;
      sqrunave[0][i] /= REDO;
      if ( bcycle!=0 ) {
	runave[1][i] /= REDO;
	sqrunave[1][i] /= REDO;
      }
      if ( acycle!=0 ) {
	runave[2][i] /= REDO;
	sqrunave[2][i] /= REDO;
      }
      mymean = Taveave/avecnt;
      stderr = sqrunave[0][i] - 2*runave[0][i]*mymean + mymean*mymean;
      if ( stderr<0 ) { stderr = 0; }
      printf("%.4lf %lf %f %lf %lf", repval[i*REPVALS], repval[i*REPVALS+1], 
	     sqrt(stderr), repval[i*REPVALS+2], repval[i*REPVALS+3]);
      if ( bcycle!=0 ) {
	stderr = sqrunave[1][i] - runave[1][i]*runave[1][i];
	if ( stderr<0 ) { stderr = 0; }
	printf(" %lf", stderr);
      }
      if ( acycle!=0 ) {
	stderr = sqrunave[2][i] - runave[2][i]*runave[2][i];
	if ( stderr<0 ) { stderr = 0; }
	printf(" %lf", stderr);
      }
      printf("\n");
    }
  }
  /*
   *   final reporting
   */
  yaps_message("Run ave T:");
  if ( DIM>20 && verbose ) 
    for (i=0; i<DIM; i++)
      yaps_message(" %f", taveave[i]/avecnt);
  yaps_message(" = %f\n", Taveave/avecnt);
  if ( bcycle!=0 )
    yaps_message("Run ave b: = %f\n", baveave/avecnt);
  if ( acycle!=0 )
    yaps_message("Run ave a: = %f\n", aaveave/avecnt);
  
  S_free(ST);
  if ( printmean ) {
    for (i=0; i<REPAVE; i++) {
      free(runave[i]);
      free(sqrunave[i]);
    }
    free(repval);
  }
  rng_free(rng);
  return 0;
}
