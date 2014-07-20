/*
 * caching lgamma + fast f(x) = lgamma(x+b)-lgamma(b)
 * Copyright (C) 2010 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 * To compile testers, use:
 *   cc -DMAINTEST -o lgamma lgamma.c -lgsl -lgslcblas -lm
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "lgamma.h"
#include "digamma.h"

/*
 *    cache all calls to lgamma()
 */

void gcache_init(struct gcache_s *lgp, double p) {
  lgp->par = p;
  lgp->lgpar = lgamma(p);
  memset(lgp->cache, 0, GCACHE*sizeof(lgp->cache[0]));
}

double gcache_value(struct gcache_s *lgp, int j) {
  if ( j<=0 )
    return 0;
  if ( j>=GCACHE )
    return lgamma(j+lgp->par) - lgp->lgpar;
  if ( lgp->cache[j]==0 ) {
    if ( j==1 )
      lgp->cache[j] = log(lgp->par);
    else if ( j==2 )
      lgp->cache[j] = log(lgp->par*(lgp->par+1));
    else if ( j==3 )
      lgp->cache[j] = log(lgp->par*(lgp->par+1)*(lgp->par+2));
   else
    lgp->cache[j] = lgamma(j+lgp->par) - lgp->lgpar;
  }
  return lgp->cache[j];
}

void pcache_init(struct gcache_s *lgp, double p) {
  lgp->par = p;
  lgp->lgpar = digamma(p);
  memset(lgp->cache, 0, GCACHE*sizeof(lgp->cache[0]));
}

double pcache_value(struct gcache_s *lgp, int j) {
  if ( j<=0 )
    return 0;
  if ( j>=GCACHE )
    return digamma(j+lgp->par) - lgp->lgpar;
  if ( lgp->cache[j]==0 ) {
    if ( j==1 )
      lgp->cache[j] = 1/lgp->par;
    else if ( j==2 ) 
      lgp->cache[j] = 1/lgp->par + 1/(1+lgp->par);
    else if ( j==3 ) 
      lgp->cache[j] = 1/lgp->par + 1/(1+lgp->par) + 1/(2+lgp->par);
    else
      lgp->cache[j] = digamma(j+lgp->par) - lgp->lgpar;
  }
  return lgp->cache[j];
}

void qcache_init(struct gcache_s *lgp, double p) {
  lgp->par = p;
  if ( p>0 ) 
    lgp->lgpar = lgamma(1-2*p) - lgamma(1-p);
  else
    lgp->lgpar = 0;
  memset(lgp->cache, 0, GCACHE*sizeof(lgp->cache[0]));
}

/*
 *  this is  S^{n+1}_{2,a}/S^n_{1,a}
 *
 *  for a>0, precompute  lga0 = lgamma(1-2*a)-lgamma(1-a)
 */
static double qval(double a, int n, double lga0) {
  if ( a<0.02 ) {
    return digamma(n+1-a) - digamma(1-a);
  }
  return (1.0 - exp(lgamma(n+1-2*a)-lgamma(n+1-a)-lga0))/a;
}
double qcache_value(struct gcache_s *lgp, int j) {
  if ( j<=0 )
    return 0;
  if ( j>=GCACHE ) 
    return qval(lgp->par,j,lgp->lgpar);
  if ( lgp->cache[j]==0 ) {
    if ( j==1 )
      lgp->cache[j] = 1/(1-lgp->par);
    else if ( j==2 ) 
      lgp->cache[j] = 3/(2-lgp->par);
    else if ( j==3 ) 
      lgp->cache[j] = (11-7*lgp->par)/(3-lgp->par)/(2-lgp->par);
    else
      lgp->cache[j] = qval(lgp->par,j,lgp->lgpar);
  }
  return lgp->cache[j];
}

/*
 *   faster lgamma(N+alpha)-lgamma(alpha)
 */
#ifndef LS_NOPOLYGAMMA
#define FDIM 2000
static double fg[FDIM];
static double fp0[FDIM];
static double fp1[FDIM];
static double fp2[FDIM];
static double fp3[FDIM];
static int fset = 0;

static void diffset(void) {
  int i;
  for (i=3; i<FDIM; i++) {
    fg[i] = lgamma(i);
    fp0[i] = digamma(i);
    fp1[i] = trigamma(i);
    fp2[i] = tetragamma(i);
    fp3[i] = pentagamma(i);
  }
  fset = 1;
}
#endif
/*
 *    faster version of lgamma(N+alpha)-lgamma(alpha)
 */
double gammadiff(int N, double alpha, double lga) {
#ifndef LS_NOPOLYGAMMA
  double val;
  if ( fset==0 ) diffset();
#endif
  if ( N<=3 ) {
    if ( N<=1 ) {
      if ( N==0 )
	return 0;
      return log(alpha);
    } 
    if ( N==2 )
      return log(alpha*(1+alpha));
    return log(alpha*(1+alpha)*(2+alpha));
  } 
#ifdef LS_NOPOLYGAMMA
  return lgamma(N+alpha) - lgamma(alpha);
#else
  else if ( alpha>0.5 ) {
    return lgamma(N+alpha) - lgamma(alpha);
  } else if ( lga != 0 ) {
    /*
     *  lgamma(alpgha) precomputed as lga
     */
    if ( N>= FDIM ) 
      val = lgamma(N+alpha);
    else
      val = fg[N] 
	+ alpha*(fp0[N]
		 + alpha/2 * (fp1[N] + alpha/3 * fp2[N]));
    return val - lga;
  }
  // lgamma(3+alpha)-lgamma(alpha)
  // + lgamma(N+alpha)
  // - lgamma(3+alpha);
  val = log(alpha*(1+alpha)*(2+alpha));
  if ( N>= FDIM ) 
    val += lgamma(N+alpha) - 
      ( fg[3] + alpha*(fp0[3] + 
		       alpha/2 * ( fp1[3] + alpha/3 * fp2[3])));
  else
    val += fg[N]-fg[3] 
      + alpha*( (fp0[N]-fp0[3]) 
		+ alpha/2 * ((fp1[N]-fp1[3])
			     + alpha/3 * (fp2[N]-fp2[3])));
  return val;
#endif
}
/*
 *    faster version of digamma(N+alpha)-digamma(alpha)
 */
double psidiff(int N, double alpha, double pa) {
#ifndef LS_NOPOLYGAMMA
  double val;
  if ( fset==0 ) diffset();
#endif
  assert(alpha>0);
  if ( N<=3 ) {
    if ( N<=1 ) {
      if ( N==0 )
	return 0;
      return 1/alpha;
    } 
    if ( N==2 )
      return 1/alpha + 1/(1+alpha);
    return 1/alpha + 1/(1+alpha) + 1/(2+alpha);
  } 
#ifdef LS_NOPOLYGAMMA
  return digamma(N+alpha) - ((pa>0)?pa:digamma(alpha));
#else
else if ( alpha>0.5 ) {
    return digamma(N+alpha) - ((pa>0)?pa:digamma(alpha));
  } else if ( pa != 0 ) {
    /*
     *  lgamma(alpgha) precomputed as pa
     */
    if ( N>= FDIM ) 
      val = digamma(N+alpha);
    else
      val = fp0[N] + alpha*(fp1[N]
			    + alpha/2 * (fp2[N] + alpha/3 * fp3[N]));
    return val - pa;
  }
  // psi(3+alpha)-psi(alpha)
  // + psi(N+alpha)
  // - psi(3+alpha);
  val = 1/alpha + 1/(1+alpha) + 1/(2+alpha);
  if ( N>= FDIM ) 
    val += lgamma(N+alpha) - 
      ( fp0[3] + alpha*(fp1[3] + 
			alpha/2 * ( fp2[3] + alpha/3 * fp3[3])));
  else
    val += fp0[N]-fp0[3] + alpha*(fp1[N]-fp1[3]) 
      + alpha*alpha/2 * (fp2[N]-fp2[3])
      + alpha*alpha*alpha/6 * (fp3[N]-fp3[3]);
  return val;
#endif
}

#ifdef MAINTEST
#define REPS 5000000
int main(int argc, char* argv[]) {
  int i;
  clock_t t1 = 0;
  double alpha = 0.2;
  double lga = lgamma(alpha);
  double pa = digamma(alpha);
  double try = 0;

  for (i=3; i<200; i+=100)
    for (alpha=0.1; alpha<1; alpha+=0.1) {
      double l1 = lgamma(i+alpha)-lgamma(alpha);
      double l2 = gammadiff(i,alpha,0);
      if ( fabs(l1-l2)/l1 > 1e-7 ) 
	printf("%d %lf %lf %lf\n", 
	       i, alpha,
	       lgamma(i+alpha)-lgamma(alpha), 
	       gammadiff(i,alpha,0) );
    }
#ifdef TIMETEST
  diffset();

  t1 = clock();
  for (i=0; i<REPS; i++) {
    try += gammadiff(100, alpha, lga);
  }
  printf("gammadiff(lga) = %.8lg\n",  
	 ((double)(clock()-t1))/CLOCKS_PER_SEC/REPS);

  t1 = clock();
  for (i=0; i<REPS; i++) {
    try += gammadiff(100, alpha, 0);
  }
  printf("gammadiff(0) = %.8lg\n",  
	 ((double)(clock()-t1))/CLOCKS_PER_SEC/REPS);

  t1 = clock();
  for (i=0; i<REPS; i++) {
    try += lgamma(100+alpha) - lga;
  }
  printf("gamma fnct = %.8lg\n",  
	 ((double)(clock()-t1))/CLOCKS_PER_SEC/REPS);

  t1 = clock();
  for (i=0; i<REPS; i++) {
    try += psidiff(100, alpha, pa);
  }
  printf("psidiff(pa) = %.8lg\n",  
	 ((double)(clock()-t1))/CLOCKS_PER_SEC/REPS);

  t1 = clock();
  for (i=0; i<REPS; i++) {
    try += psidiff(100, alpha, 0);
  }
  printf("psidiff(0) = %.8lg\n",  
	 ((double)(clock()-t1))/CLOCKS_PER_SEC/REPS);

  t1 = clock();
  for (i=0; i<REPS; i++) {
    try += digamma(100+alpha) - pa;
  }
  printf("psi fnct = %.8lg\n",  
	 ((double)(clock()-t1))/CLOCKS_PER_SEC/REPS);
#endif

  return 0;
}
#endif
