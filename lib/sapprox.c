/*
 * Exact S computation
 * Copyright (C) 2009 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *  Unoptimised, no attempt to cache any values.
 *
 *  Exact calculations for S, for theory, see:
 *      http://arxiv.org/abs/1007.0296
 *          
 */
#include <assert.h>
#include <math.h>
#include "digamma.h"

/*
 *    approx calculation for M<=4,
 *    but is exact when a==0, and OK estimate
 *    for a small
 */
double S_approx(int n, int m, float a) {
  assert(a>0);
  if ( n==m ) 
    return 0.0;
  if ( n<m ) 
    return -HUGE_VAL;
  if ( m==1 ) {
    return lgamma(n-a) - lgamma(1-a);
#ifndef LS_NOPOLYGAMMA
  } else if ( a<0.001 ) {
    if (m==2 ) {
      return lgamma(n-a) - lgamma(1-a) + 
	log(digamma(n-a) - digamma(1-a));
    } else if ( m==3 ) {
      double gg = digamma(n-a) - digamma(1-a);
      return lgamma(n-a) - lgamma(1-a) - log(2.0) + 
	log(trigamma(n-a) - trigamma(1-a) + gg*gg);
    } else if ( m==4 ) {
      double gg = digamma(n-a) - digamma(1-a);
      return lgamma(n-a) - lgamma(1-a) - log(6.0) + 
	log((tetragamma(n-a) - tetragamma(1-a)) +
	    3*(trigamma(n-a) - trigamma(1-a))*gg + gg*gg*gg);
    }
#endif
  } else if ( m==2 ) {
    double ga = lgamma(n-a) - lgamma(1-a);
    double g2a = lgamma(n-2*a) - lgamma(1-2*a);
    return g2a - log(a) + log(exp(ga-g2a)-1.0);
  } else if ( m==3 ) {
    double ga = lgamma(n-a) - lgamma(1-a);
    double g2a = lgamma(n-2*a) - lgamma(1-2*a);
    double g3a = lgamma(n-3*a) - lgamma(1-3*a);
    return g3a - 2*log(a) - log(2.0) +
      log(exp(ga-g3a)-2*exp(g2a-g3a)+1.0);
  } else if ( m==4 ) {
    double ga = lgamma(n-a) - lgamma(1-a);
    double g2a = lgamma(n-2*a) - lgamma(1-2*a);
    double g3a = lgamma(n-3*a) - lgamma(1-3*a);
    double g4a = lgamma(n-4*a) - lgamma(1-4*a);
    return g4a - 3*log(a) - log(6.0) +
      log(exp(ga-g4a)-3*exp(g2a-g4a)+3*exp(g3a-g4a)-1.0);
  }
  return -HUGE_VAL;
}

/*
 *    approx derivative calculation for M<=4
 */
double S_approx_da(int n, int m, float a) {
  double snm;
  assert(a>0);
  if ( n==m ) 
    return 0.0;
  if ( n<m ) 
    return -HUGE_VAL;
  if ( m==1 ) {
    return -(digamma(n-a) - digamma(1-a));
  } 
  snm = S_approx(n, m, a);
  if ( m==2 ) {
    double ga = lgamma(n-a) - lgamma(1-a);
    double g2a = lgamma(n-2*a) - lgamma(1-2*a);
    double dga = -(digamma(n-a) - digamma(1-a));
    double dg2a = -2.0*(digamma(n-2*a) - digamma(1-2*a));
    return (exp(ga-snm)*dga - exp(g2a-snm)*dg2a - 1)/a;
  } else if ( m==3 ) {
    double ga = lgamma(n-a) - lgamma(1-a);
    double g2a = lgamma(n-2*a) - lgamma(1-2*a);
    double g3a = lgamma(n-3*a) - lgamma(1-3*a);
    double dga = -(digamma(n-a) - digamma(1-a));
    double dg2a = -2*(digamma(n-2*a) - digamma(1-2*a));
    double dg3a = -3*(digamma(n-3*a) - digamma(1-3*a));
    return  - 2/a + (exp(ga-snm)*dga-2*exp(g2a-snm)*dg2a+exp(g3a-snm)*dg3a)/2/a/a;
  } else if ( m==4 ) {
    double ga = lgamma(n-a) - lgamma(1-a);
    double g2a = lgamma(n-2*a) - lgamma(1-2*a);
    double g3a = lgamma(n-3*a) - lgamma(1-3*a);
    double g4a = lgamma(n-4*a) - lgamma(1-4*a);
    double dga = -(digamma(n-a) - digamma(1-a));
    double dg2a = -2*(digamma(n-2*a) - digamma(1-2*a));
    double dg3a = -3*(digamma(n-3*a) - digamma(1-3*a));
    double dg4a = -4*(digamma(n-4*a) - digamma(1-4*a));
    return - 3/a + 
      (exp(ga-snm)*dga-3*exp(g2a-snm)*dg2a+3*exp(g3a-snm)*dg3a-exp(g4a-snm)*dg4a)/3/a/a/a;
  }
  return -HUGE_VAL;
}

