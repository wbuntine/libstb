/*
 * Sampling hyperparameters
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

#ifndef __PSAMPLE_H
#define __PSAMPLE_H

#include "stable.h"
#include "srng.h"
#include "lgamma.h"

/*
 *  Two methods of sampling a
 *      * rebuilding S table with each change of a
 *      * sample table totals to get S-free likelihood
 *  defining this does the latter, usually 70% faster
 *  but its an experimental approach and may have trouble if
 *  table counts become very large
 */
// #define SAMPLEA_M
/*
 *     Define to use ARS rather than slice sampler
 *     NOTE:  if you have disabled polygamma in "digamma.h"
 *            then you must use ARS since the slice
 *            sampler requires the use of digammaInv()
 */
#define PSAMPLE_ARS
#ifdef PSAMPLE_ARS
#include "arms.h"
#else
/*
 *    return non-zero on fatal error;
 *    posterior must be unimodal since it fiddles bounds
 */
int SliceSimple(double *xp,      /*  input and output result */
		double (*post)(double, void *), 
		double *bounds,  /*  left and right bound */
		rngp_t rng, 
		int loops,   /* loops of slice sampler before returning */ 
		void *pars   /* args for post()  */
		);
#endif

/*
 *    bounds on concentration, b
 *    NB.   slight problem, usually b>-a !!
 */
#define B_MIN 0.01
#define B_MAX 2000

/*
 *  type used for larger counts
 */
typedef uint32_t scnt_int;
/*
 *  type used for table counts (maybe smaller)
 */
typedef uint16_t stcnt_int;

/*
 *   b_in : its an MCMC sampler, so input value needed
 *   I : dimension of N[] and T[]
 *   shape, scale : gamma prior parameters for b
 *   N[], T[] : totals (over topics) for document counts and tables
 *   apar : discount a
 *   loops : of the slice sampler before returning
 *   verbose : print more details if >1
 */
double sampleb(double b_in, 
	       int I,
	       double shape, double scale, 
	       scnt_int *N, scnt_int *T, 
	       double apar, 
	       rngp_t rng, int loops, int verbose);

/*
 *    bounds on discount, a
 */
#define A_MIN 0.01
#define A_MAX 0.98
/*
 *    don't allow shifts bigger than this in any one MCMC step
 */
#define SQUEEZEA 0.2

/*
 *  n, t counts have IxK[i] dimension (second dimension varies with i=0,..,K-1
 *  T[i] is the total over k=0,..,K[i] for t[i][k]
 *  bpar[i] gives b value to use for each i=0,..,K-1
 *  actual n,t values can be in matrices passed, or can be
 *      obtained using the callback getval(&n,&t,i,k)
 *  Creates its own S table internally
 */
double samplea(double apar, 
	       int I, int *K,     /* dimensions */
	       scnt_int *T,           /* topicwise totals of t[][] */
	       scnt_int **n, stcnt_int **t, /* counts and tables */
	       void (*getval)(scnt_int *n, stcnt_int *t, unsigned i, unsigned k),
	       double *bpar,
	       rngp_t rng, int loops, int verbose);

double samplea2(double mya, stable_t  *S,
		int I, int *K, scnt_int *T,
		scnt_int **n, stcnt_int **t,
		void (*getval)(scnt_int *n, stcnt_int *t, unsigned i, unsigned k),
		double *bpar,
		rngp_t rng, int loops, int verbose);
#endif
