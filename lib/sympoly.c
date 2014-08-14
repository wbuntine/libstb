/*
 * Symmetric polynomials
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
 *  See http://en.wikipedia.org/wiki/Newton%27s_identities
 *  for math details
 */

#include <stdio.h>  
#include <unistd.h>   
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

/*
 *   change this to use another RNG
 */
#include "srng.h"
#include "sympoly.h"

/*
`*  Compute elementary symmetric polynomials
 *     res[0] = 1
 *     res[1] = e_1(x_0,...,x_{K-1})
 *     res[2] = e_2(x_0,...,x_{K-1})
 *     ...
 *     res[K] = e_K(...)
 *
 *    use an alternative recursion:
 *      e_H(x_0,...,x_{K-1}) = 
 *            1_{H<K} && e_H(x_0,...,x_{K-2}) + x_{K-1}e_{H-1}(x_0,...,x_{K-2})
 *    build values in situ in the res[]
`*
 *    but convert this to better behaved form
 *
 *      f_H(x_0,...,x_{K-1}) == \prod_{k=0}{^{K-1} (1/x_k)^{x_k>1}
 *                               * e_H(x_0,...,x_{K-1})
 *    which has recursion:
 *      f_H(x_0,...,x_{K-1}) = 
 *         1_{H<K} && f_H(x_0,...,x_{K-2}) + x_{K-1}f_{H-1}(x_0,...,x_{K-2})
 *                               if x_{K-1}<=1
 *         1_{H<K} && f_H(x_0,...,x_{K-2})/x_{K-1} + f_{H-1}(x_0,...,x_{K-2})
 *                               otherwise
 *
 *    if values get too big, log overflow will appear in (*overflow),
 *    reconstruct values as
 *       res[0], exp(*overflow)*(res[1],...,res[K])
 *
 *  return non-zero on error
 */
int sympoly(int K, int BK, double *val, double *res, double *overflow) {
  int h, k;
  double mult = 1;           //  keeps scaled res[0]
  /*
   *   e_0=1 regardless
   */
  BK++;
  if ( BK>K )
    BK = K;
  if ( BK<=1 )
    BK = 2;
  *overflow = 0;
  res[0] = 1;
  for (k=1; k<=K; k++) 
    res[k] = 0;
  if ( K==0 )
    return 0;
  res[1] = val[0];
  for (k=1; k<K; k++) {
    double vk = val[k];
    if ( vk>1 ) {
      double lastv = mult;
      mult /= vk;
      *overflow += log(vk);
      for (h=1; h<=k && h<BK; h++) {
	double thisv = res[h];
	res[h] = thisv/vk + lastv;
	lastv = thisv;
      }
      res[h] = lastv;
    } else {
      double lastv = mult;
      for (h=1; h<=k && h<BK; h++) {
	double thisv = res[h];
	res[h] += vk*lastv;
	lastv = thisv;
      }
      res[h] = vk*lastv;
    }
  }
  if ( *overflow<15 ) {
    /*
     *  remove overflow if doesn't appear useful
     */
    double eo = exp(*overflow);
    for (h=1; h<=BK; h++)
      res[h] *= eo;
    *overflow = 0;
  }
  return 0;
}

/*
 *    sample exactly H of the K features proportional to
 *    occurrences in the lem.sym.poly.
 *
 *    use an alternative recursion:
 *      e_H(x_0,...,x_{K-1}) = 
 *            1_{H<K} && e_H(x_0,...,x_{K-2}) + x_{K-1}e_{H-1}(x_0,...,x_{K-2})
 *    build all the values in a table, and then recurse back on K using:
 *
 *    if H==K, then all remaining on
 *    if H==0, then none remaining on
 *    else
 *       p(x_{K} off) = e_{H}(x_0,...,x_{K-1})/e_H(x_0,...,x_{K})
 *       if x_{K} on then decrease H
 *
 *    return 0 if OK, non-zero on error
 */
static uint32_t sympoly_sample_full(int K, int H, double *val, rngp_t rng) {
  int h, k;
  /*
   *    esp[k+(h-1)*K] =  f_{h}(x_0,...,x_k)  for h>0
   *    i.e.,   don't tabulate h=0 case, is unused
   */
  double *esp, espp[SYMPOLY_MAX*SYMPOLY_MAX];
  uint32_t wr;
  double mult = 1;
  assert(K>1);
  assert(H>1);
  assert(H<K);
  if ( H*K<SYMPOLY_MAX*SYMPOLY_MAX )
    esp = &espp[0];
  else 
    esp = malloc(sizeof(*esp)*K*H);
  if ( !esp )
    return 0;
  /*
   *   fill up matrix, O(HK), 
   *   cannot do the vector fill as before since we need
   *   all values for all values of H and K
   */
  esp[0] = val[0];
  for (k=1; k<K; k++) {
    double vk = val[k];
    if ( vk>1 ) {
      double lastv = mult;
      mult /= vk;
      for (h=0; h<k && h<H; h++) {
	double thisv = esp[k-1+h*K];
	esp[k+h*K] = thisv/vk + lastv;
        lastv = thisv;
      }
      if ( h==k )
	esp[k+h*K] = lastv;
    } else {
      double lastv = mult;
      for (h=0; h<k && h<H; h++) {
	double thisv = esp[k-1+h*K];
	esp[k+h*K] = thisv + vk*lastv;
        lastv = thisv;
      }
      if ( h==k )
	esp[k+h*K] = vk*lastv;
    }
  }
#if 0
  printf("Comp(k=%d): ", K) ;
  for (h=0; h<H; h++)
    printf(" %lg", esp[K-1+h*K]);
  printf("\n");
#endif

  /*
   *    recurse back
   */
  wr = 0;
  for (k=K-1, h=H;  k>=h && h>0;  k--) {
    assert(k>=1);
    if ( val[k]<=1 ) {
      if ( esp[k+(h-1)*K]* rng_unit(rng) >= esp[k-1+(h-1)*K] ) {
	wr |= 1U<<(unsigned)k;
	h--;
      }
    } else {
      if ( esp[k+(h-1)*K]* rng_unit(rng)  >=  esp[k-1+(h-1)*K]/val[k] ) {
	wr |= 1U<<(unsigned)k;
	h--;
      }
    }
  }
  if ( esp!=&espp[0] )
    free(esp);
  assert(h==(k+1) || h==0);
  if ( h>0 )
    wr |= (1U<<(unsigned)h) - 1;
  return wr;
}
/*

 *    sample exactly H of the K features proportional to
 *    occurrences in the elem.sym.poly.
 *
 *    return 0 if OK, non-zero on error
 */
uint32_t sympoly_sample(int K, int H, double *val, rngp_t rng) {
  int k ;
  double sum;
  if ( H>K || K==0 || H==0 )
    return 0U;
  if ( H==K ) {
    return (1U<<(unsigned)H) - 1;
  }
  if ( H>1 ) { 
    return sympoly_sample_full(K, H, val, rng);
  }
  assert(H==1);
  /*
   *    standard vector sampler for H=1 case
   */
  sum = 0;
  for (k=0; k<K; k++) {
    sum += val[k];
  }
  sum *= rng_unit(rng);
  for (k=0; k<K && sum>0; k++)
    sum -= val[k];
  assert(k>0);
  return 1U<<(unsigned)(k-1U);
}

/*
 *   revision 49 in SVN has the old version concurrent
 *   which was used to test if this worked OK
 */
//#define TEST
#ifdef TEST
#include <time.h>

#define KMAX 20
int main(int argc, char* argv[]) {
  int K=10;
  double val[KMAX+1], 
    res[KMAX+1];
  double overflow;
  int k;
  rngp_t rng=NULL;

  // srand48(time(NULL));
  srand48(10);

  for (k=0; k<K; k++) val[k] = drand48();
  for (k=0; k<K; k++) if ( drand48()>0.5 ) val[k] = 1.0/val[k];

  printf("Val: ");
  for (k=0; k<K; k++)  printf(" %lg", val[k]);
  printf("\n");

#if 0
  sympoly_old(K, val, res);
  printf("Res : ");
  for (k=0; k<=K; k++)  printf(" %lg", res[k]);
  printf("\n");
#endif
  sympoly(K, K, val, res, &overflow);
  printf("Res overflow = %lg\n", overflow);
  printf("Res2:  %lg", res[0]);
  for (k=1; k<=K; k++)  printf(" %lg", exp(overflow)*res[k]);
  printf("\n");
  printf("Res2:  %lg", res[0]);
  for (k=1; k<=K; k++)  printf(" %lg", res[k]);
  printf("\n");
  sympoly(K, K/2, val, res, &overflow);
  printf("Res bound overflow = %lg\n", overflow);
  printf("Res bound:  %lg", res[0]);
  for (k=1; k<=K/2; k++)  printf(" %lg", exp(overflow)*res[k]);
  printf("\n");

  
  printf("Sampled2: ");
  for (k=0; k<10; k++) 
    printf(" %x",  sympoly_sample(K, 6, val, rng));
  printf("\n");
  return 1;
}
#endif
