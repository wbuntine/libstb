/*
 *    A simple interface to the chosen RNG
 *
 *    currently using rand48() so it has no local
 *    data structures, moreover it *cannot* be used in
 *    multi-threaded programs!!!
 */
#ifndef __RNG_H
#define __RNG_H

#include <time.h>

/*
 *  these are defined in  gslrandist.c;
 *  replace with your own if needed
 */
double gsl_rng_gaussian_ziggurat (const double sigma);
double gsl_rng_beta (const double a, const double b);
double gsl_rng_gamma (const double a);

typedef void *rngp_t;

/*
 *  these macros must be provided;
 *  for GSL or similar, the rng variable would be the local
 *  strucure
 */
#define rng_seed(rng,seed) srand48(seed);
#define rng_time(rng,seed)  {  *(seed) = time(NULL); srand48(*(seed)); }
#define rng_unit(rng) drand48()
#define rng_beta(rng,a,b) gsl_rng_beta(a,b)
#define rng_gamma(rng,a) gsl_rng_gamma(a)
#define rng_gaussian(rng,a) gsl_rng_gaussian_ziggurat(a)
#define rng_free(rng)  
#endif
