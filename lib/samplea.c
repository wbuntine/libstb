/*
 * Sampling discount parameter for Pitman-Yor using Slice sampler
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
 *    This is rather naive.   We should probably sample the
 *    table assignments, i.e., figure out a partition of each
 *    n across its t tables, and then the posterior is just
 *    gamma functions, not Stirling numbers.
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "srng.h"
#include "psample.h"
#include "stable.h"

typedef struct aldata {
  int    I;
  int    *K;
  scnt_int    **n;
  stcnt_int    **t;
  void (*val)(scnt_int *n, stcnt_int *t, unsigned i, unsigned k);
  int    maxt;
  int    maxn;
  scnt_int    *T;
#ifdef SAMPLEA_M
  stcnt_int *m;
#endif
  double    *bpar;
  stable_t  *S;
  int verbose;
} ALData;

static double aterms(double x, void *mydata) {
  int i, k;
  ALData *mp = (ALData *)mydata;
  double val = 0;
  if ( x<=0 ) {
    fprintf(stderr,"Illegal discount value in aterms()\n");
    exit(1);
  }
  if ( mp->verbose>1 ) {
    fprintf(stderr,"Extending S for M=%d a=%lf\n", mp->maxt,x);
  }
  if ( mp->S )
    S_remake(mp->S,x);
  else
    mp->S = S_make(mp->maxn,mp->maxt,mp->maxn,mp->maxt,x,S_STABLE);
  if ( !mp->S ) {
    fprintf(stderr,"Out of memory for S table\n");
    exit(1);
  }
  for (i=0; i<mp->I; i++) {
    val += mp->T[i] * log(x) + 
      lgamma(mp->T[i]+mp->bpar[i]/x) - lgamma(mp->bpar[i]/x);
    if ( mp->val ) {
      for (k=0; k<mp->K[i]; k++) {
	scnt_int n;
	stcnt_int t;
	mp->val(&n, &t, i, k);
	if ( n>1 )
	  val += S_S(mp->S, n, t);
      }
    } else {
      for (k=0; k<mp->K[i]; k++)
	if ( mp->n[i][k]>1 )
	  val += S_S(mp->S, mp->n[i][k],mp->t[i][k]);
    }
  }
  return val;
}

#ifdef SAMPLEA_M
#define LGCACHE   
static double aterms2(double x, void *mydata) {
  int i, k;
  ALData *mp = (ALData *)mydata;
  double val = 0;
#ifdef LGCACHE
  struct gcache_s lgp;
#else
  double lg1a = lgamma(1-x);
#endif
  stcnt_int *mm = mp->m;
  if ( x<=0 ) {
    fprintf(stderr,"Illegal discount value in aterms2()\n");
    exit(1);
  }
#ifdef LGCACHE
  gcache_init(&lgp, 1-x);
#endif
  for (i=0; i<mp->I; i++) {
    val += mp->T[i] * log(x) + 
      lgamma(mp->T[i]+mp->bpar[i]/x) - lgamma(mp->bpar[i]/x);
    for (k=0; k<mp->K[i]; k++) {
      scnt_int n;
      stcnt_int t;
      if ( mp->val ) {
	mp->val(&n, &t, i, k);
      } else {
	n = mp->n[i][k];
	t = mp->t[i][k];
      }
      if ( n>0 ) {
	if ( t==n ) {
	  ;
	} else if ( t==1 ) {
#ifdef LGCACHE
	  val += gcache_value(&lgp, n-1);
#else
	  val += lgamma(n-x) - lg1a; 
#endif
	} else {
	  int l;
	  for (l=t-2; l>=0; l--) {
	    if ( mm[l]>1 ) 
#ifdef LGCACHE
	      val += gcache_value(&lgp, mm[l]-1);
#else
	      val += lgamma(mm[l]-x) - lg1a; 
#endif
	    n -= mm[l];
	  }
	  assert(n>0);
	  if ( n>0 ) 
#ifdef LGCACHE
	    val += gcache_value(&lgp, n-1);
#else
	    val += lgamma(n-x) - lg1a; 
#endif
	  mm += t-1;
	}
      }
    }
  }
  return val;
}
#endif

/*
 *    use prior a is uniform
 */
double samplea(double mya, 
	       int I, int *K, scnt_int *T,
	       scnt_int **n, stcnt_int **t,
	       void (*getval)(scnt_int *n, stcnt_int *t, unsigned i, unsigned k),
	       double *bpar,
	       rngp_t rng, int loops, int verbose) {
  double inita[3] = {A_MIN,1,A_MAX};
  int i, k;
  ALData ald;
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
  ald.T = T;
  ald.n = n; 
  ald.t = t;
  ald.I = I;
  ald.K = K;
  ald.val = getval;
  ald.bpar = bpar;
  ald.verbose = verbose;
  ald.maxt = 1; 
  ald.maxn = 1; 
  ald.S = NULL;
  if ( getval ) {
    for (i=0; i<I; i++)
      for (k=0; k<K[i]; k++) {
	scnt_int myn;
	stcnt_int myt;
	getval(&myn, &myt, i, k);
	if ( myt>=ald.maxt )
	  ald.maxt = myt+1;
	if ( myn>=ald.maxn )
	  ald.maxn = myn+1;
      }
  } else {
    for (i=0; i<I; i++)
      for (k=0; k<K[i]; k++) {
	if ( t[i][k]>=ald.maxt )
	  ald.maxt = t[i][k]+1;
	if ( n[i][k]>=ald.maxn )
	  ald.maxn = n[i][k]+1;
      }
  }
#ifdef PSAMPLE_ARS
  arms_simple(3,  inita, inita+2, aterms,
	      &ald, 0, inita+1, &mya);
  if ( mya<inita[0] ||  mya>inita[2] ) {
    fprintf(stderr,"Arms_simple(apar) returned value out of bounds\n");
    exit(1);
  }
#else
  inita[1] = A_MAX;
  if ( SliceSimple(&mya, aterms, inita, rng, loops, &ald) ) {
    fprintf(stderr,"SliceSimple error\n");
    exit(1);
  }
#endif
  if ( ald.S ) S_free(ald.S);
  return mya;
}

#ifdef SAMPLEA_M

/*
 *  log(exp(x)-exp(y)) =
 *     x + log(1-exp(y-x))
 */
double logminus(double x, double y) {
  if ( y>=x ) 
    return -HUGE_VAL;
  if ( y-x<-80 ) 
    return x - exp(y-x);
  return x + log(1-exp(y-x));
}

/*
 *    use prior a is uniform
 */
double samplea2(double mya, stable_t  *S,
		int I, int *K, scnt_int *T,
		scnt_int **n, stcnt_int **t,
		void (*getval)(scnt_int *n, stcnt_int *t, unsigned i, unsigned k),
	       double *bpar,
	       rngp_t rng, int loops, int verbose) {
  double inita[3] = {A_MIN,1,A_MAX};
  int i, k;
  ALData ald;
  stcnt_int *mp;
  int n_m = 0;
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
  ald.T = T;
  ald.n = n; 
  ald.t = t;
  ald.I = I;
  ald.K = K;
  ald.val = getval;
  ald.bpar = bpar;
  ald.verbose = verbose;
  ald.maxt = 1; 
  ald.maxn = 1; 
  ald.S = NULL;
  n_m = 0;
  for (i=0; i<I; i++)
    for (k=0; k<K[i]; k++) {
      if ( t[i][k]>1 && t[i][k]<n[i][k]) 
	n_m += t[i][k]-1 ;
    }
  ald.m = malloc(sizeof(*ald.m)*n_m);
  if ( !ald.m ) {
    fprintf(stderr,"Out of memory for samplea()\n");
    exit(1);
  }
  // fprintf(stderr,"m[%d][%d] tabl space = %d\n", I, K[0], n_m);
  mp = ald.m;
  for (i=0; i<I; i++)
    for (k=0; k<K[i]; k++) {
      if ( t[i][k]>1 && t[i][k]<n[i][k]) {
	int N = n[i][k];
	double ptot = S_S(S,N,t[i][k]);
	double rem = ptot + log(rng_unit(rng));
	int M;
	// fprintf(stderr,"Sampling[%d][%d] n=%d, t=%d: ", i,k, N, t[i][k]);
	for (M=t[i][k]-1; M>=1; M--) {
	  //  each round instantiates another count
	  int l;
	  double fact = 0.0;
	  for (l=1; l<=N-M; l++) {
	    double term;
	    if ( l>1 ) 
	      fact += log((l-mya)*(N-l+1)/(l-1));
	    term = fact+S_S(S,N-l,M)-ptot;
	    if ( term>= rem )
	      break;
	    rem = logminus(rem, term);
	  }
	  if ( l>N-M ) l = N-M;
	  mp[M-1] = l;
	  N -= l;
	  assert(N>=1);
	  // fprintf(stderr," %d", l);
	}
	// fprintf(stderr," %d\n", N);
	mp += t[i][k]-1;
      } 
    }
#ifdef PSAMPLE_ARS
  arms_simple(3,  inita, inita+2, aterms2,
	      NULL, 0, inita+1, &mya);
  if ( mya<inita[0] ||  mya>inita[2] ) {
    fprintf(stderr,"Arms_simple(apar) returned value out of bounds\n");
    exit(1);
  }
#else
  inita[1] = A_MAX;
  if ( SliceSimple(&mya, aterms2, inita, rng, loops, &ald) ) {
    fprintf(stderr,"SliceSimple error\n");
    exit(1);
  }
#endif
  if ( ald.S ) S_free(ald.S);
  return mya;
}
#endif
