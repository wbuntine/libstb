/*
 * Elementary symmetric polynomials
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
 *  for math details of elementary symmetric polynomials,
 *  but we don't use their recursions because they are unstable.
 *  
 *  Based on an alternative recursion (on K):
 *      e_H(x_0,...,x_{K-1}) = 
 *            1_{H<K} && e_H(x_0,...,x_{K-2}) + x_{K-1}e_{H-1}(x_0,...,x_{K-2})
 *  to build values in situ in a vector or matrix.
 *  this recursion is more stable and allows a simple trick to
 *  deal with overflow.  If many x_k values >1, can get overflow.
 *  So instead we collect these together and run the recursion on:
 *  f_H(x_0,...,x_{K-1}) == \prod_{k=0}{^{K-1} (1/x_k)^{x_k>1}
 *                               * e_H(x_0,...,x_{K-1})
 */

#ifndef __SYMPOLY_H
#define __SYMPOLY_H

#include <stdint.h>
#include "srng.h"

/*
 *   for optional statically declared arrays, but these are just
 *   buffers that are worked around using malloc(), so no need
 *   to change ....
 */
#define SYMPOLY_MAX 10

/*
`*  Compute elementary symmetric polynomials
 *     res[0] = 1
 *     res[1] = e_1(x_0,...,x_{K-1})
 *     res[2] = e_2(x_0,...,x_{K-1})
 *     ...
 *     res[BK] = e_BK(...)
 *  i.e.  BK used as bound, should be <=K
 *
 *  is O(K^2) when BK=K, otherwise O(K*BK)
 *
 *  if values get too big, log overflow will appear in (*overflow),
 *    reconstruct values as
 *       res[0], exp(*overflow)*(res[1],...,res[BK])
 *  otherwise  overflow=0
 *
 *  return non-zero on error (e.g.,  K too big)
 */
int sympoly(int K, int BK, double *val, double *res, double *overflow);

/*
 *    sample exactly H of the K features proportional to
 *    occurrences in the elem.sym.poly.
 *
 *    is O(HK) for H<=K
 *
 *    results returned as a bit-vector, 
 *           i.e.,  max. for H is 31!
 *
 *    return 0 if OK, non-zero on error
 */
uint32_t sympoly_sample(int K, int H, double *val, rngp_t rng);
 
#endif
