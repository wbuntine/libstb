/*
 * Sampling concentration parameter for Pitman-Yor using Slice sampler
 * Copyright (C) 2012 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "srng.h"
#include "digamma.h"
#include "psample.h"

typedef struct bldata {
  double shape;
  double Q;
  int    I;
  scnt_int    *T;
  double    apar;
} BLData;

static double bterms(double x, void *mydata) {
  BLData *mp = (BLData *)mydata;
  int i;
  double lg = lgamma(x/mp->apar);
  double val = -mp->Q * x + (mp->shape-1)*log(x);
  for (i=0; i<mp->I; i++)
    val += lgamma(mp->T[i]+x/mp->apar) - lg;
  return val;
}

#ifndef PSAMPLE_ARS
/*
 *   this is a maximiser for bterms();
 *   but we set the termination condition very weakly,
 *   we just need a high probability starting point
 */
#define B_ERROR 1.0e-4      // relative error
#define B_LOOPS 5           // maximum loops
static double bmax(double x, BLData *mp) {
  double x_prime=x;
  int loops=B_LOOPS;
  int i;
  if ( x<=0 ) {
    fprintf(stderr,"Illegal concentration value in bmax()\n");
    exit(1);
  }
  x *= 1.1;
  while ( fabs((x-x_prime)/x)> B_ERROR && --loops>0 ) {
    double val = (mp->shape-1)*mp->apar/x - mp->Q*mp->apar;
    for (i=0; i<mp->I; i++) 
      val += digamma(mp->T[i]+x/mp->apar);
    x = x_prime;
    x_prime = mp->apar * digammaInv(val/mp->I);
  }
  return x_prime;
}
#endif

/*
 *    use  prior b  is Gamma(shape,scale)
 *    auxillary var  q_i  is Beta(b,N_i)
 *    When a==0
 *        then  b|q is Gamma(\sum_i T_i+shape, 1/(1/scale+\sum_i log(1/q_i)))
 *    When a>0
 *        do sampling on posterior
 */
double sampleb(double b_in, int I, double shape, double scale, 
	       scnt_int *N, scnt_int *T, double apar,
	       rngp_t rng, int loops, int verbose) {
  double Q;
  double q;
  double myb;
  int i;
  if ( scale<=0 ) {
    fprintf(stderr,"Illegal scale in sampleb()\n");
    exit(1);
  }
  Q = 1.0/scale;
  for (i=0; i<I; i++) {
    if ( N[i]<=0 ) 
      continue;
    q = rng_beta(rng, b_in, (int)N[i]);
    if ( q<=0 ) {
      fprintf(stderr,"Illegal q in sampleb(b=%lf)\n", b_in);
      exit(1);
    }
    Q -= log(q);
  }
  if ( apar==0 ) {
    double Tsum = shape;
    for (i=0; i<I; i++) 
      Tsum += T[i];
    if ( Tsum>400 ) {
      //   gamma is near enough a Gaussian with tiny std dev.
      do {
	myb = Tsum + rng_gaussian(rng, 1) * sqrt(Tsum);
      } while ( myb<=0 );
    } else
      myb = rng_gamma(rng, Tsum);
    myb /= Q;
    if ( myb<B_MIN )
      myb = B_MIN;
    if ( myb>B_MAX )
      myb = B_MAX;
    if ( verbose>1 )
      fprintf(stderr,"Sample b ~ gamma(%lg,%lg) = %lf\n", Tsum, Q, myb);
  } else {
    double initb[3] = {B_MIN,1,B_MAX};
    BLData bld;
    bld.Q = Q;
    bld.I = I;
    bld.T = T;
    bld.apar = apar;
    bld.shape = shape;
#ifdef PSAMPLE_ARS
    initb[1] = b_in;
    if ( fabs(initb[1]-B_MAX)/B_MAX<0.00001 ) {
      initb[1] = B_MAX*0.999 + B_MIN*0.001;
    }
    if ( fabs(initb[1]-B_MIN)/B_MIN<0.00001 ) {
      initb[1] = B_MIN*0.999 + B_MAX*0.001;
    }
    arms_simple(3,  initb, initb+2, bterms, &bld, 0, initb+1, &myb);
    if ( myb<B_MIN ||  myb>B_MAX ) {
      fprintf(stderr,"Arms_simple(bpar) returned value out of bounds\n");
      exit(1);
    }
#else
    /*
     *   call the maximiser, but terminate early;
     *   to get a reasonable starting point,
     *   otherwise the slice sampler behaves poorly
     */
    myb = bmax(b_in, &bld);
    if ( verbose>1 )
      fprintf(stderr,"Max b (%lg,%lg) -> %lg\n", b_in, Q, myb);
    initb[1] = B_MAX;
    if ( SliceSimple(&myb, bterms, initb, rng, loops, &bld) ) {
      fprintf(stderr,"SliceSimple error\n");
      exit(1);
    }
#endif
    if ( verbose>1 )
      fprintf(stderr,"Sample b ~ G(%lg) = %lf\n", Q, myb);
  }
  return myb;
}





