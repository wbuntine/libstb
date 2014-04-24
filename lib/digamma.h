/*
 *  Define digamma() and polygamma() here.
 *
 *  We use the GSL special functions, GSL is available
 *  on most platforms including MAC and Windows.
 *
 *  Another option is the Cephes library with psi() and zeta()
 *           http://www.netlib.org/cephes/
 *  which have nice, small self contained functions without all the
 *  GSL cruft.
 */

#ifndef __DIGAMMA_H
#define __DIGAMMA_H

/*
 *  Define this to switch off use of polygamma.
 *  Polygamma (trigamma, etc.) is used in some
 *  functions to speed things up.  So disabling
 *  it wont stop the functions working:
 *            gammadiff(), psidiff(), S_approx()
 *  But it also disables digammaInv() which
 *  uses a Newton-Raphson step.
 */
#define LS_NOPOLYGAMMA
/*
 *  Leaving this undefined means we use the library functions grabbed
 *  from the Mathlib and under GPL.  This makes the code self
 *  contained since no other libraries will be needed.
 *  Define this if you have GSL available to use instead.
 *  In which case, you will also need to modify the Makefile
 *  to include GSL.
 */
// #define GSL_POLYGAMMA

#ifdef LS_NOPOLYGAMMA
/*
 *   Radford Neal's implementation in digamma.c
 */
double digammaRN(double x);
#define digamma(x) digammaRN(x)

#else
double digammaInv(double x);

#ifndef GSL_POLYGAMMA
/*
 *   Mathlib implementation
 */
double MLpsigamma(double x, double deriv);
double MLdigamma(double x);
double MLtrigamma(double x);
double MLtetragamma(double x);
double MLpentagamma(double x);

#define digamma(x) MLdigamma(x)
#define trigamma(x) MLtrigamma(x)
#define tetragamma(x) MLtetragamma(x)
#define pentagamma(x) MLpentagamma(x)
#else
/*
 *   GSL library
 */
#include <gsl/gsl_sf.h>

#define digamma(x) gsl_sf_psi(x)
#define trigamma(x) gsl_sf_psi_n(1,x)
#define tetragamma(x) gsl_sf_psi_n(2,x)
#define pentagamma(x) gsl_sf_psi_n(3,x)

#endif
#endif

#endif

