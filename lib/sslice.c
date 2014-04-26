/*
 * Simple slice sampler on log concave distribution
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
 *  Simple slice sampler on log-concave distribution
 *          
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "srng.h"
// #include "yaps.h"

#define TOOMANY 200
/*
 *    on input, *xp should be an approximate maxima
 *               to make sampling behave beter
 *    post() is a log posterior calc, so handle accordingly
 *    do a few loops to let burn in
 *
 *   returns 1 on error
 */
int SliceSimple(double *xp, 
		double (*post)(double, void *), 
		double *bounds,
		rngp_t rng, int loops, void *pars) {
  double x = *xp, y;
  double range[2];
  int tries;
  if ( x<bounds[0] || x>bounds[1] ) {
    fprintf(stderr,
	    "SliceSimple: input value %lf outside bounds [%lg,%lg]\n",
	    x, bounds[0], bounds[1]);
    return 1;
  }
  //  yaps_message("SliceSimple: x=%lg\n", x);
  while ( loops-->0 ) {
    y = post(x,pars);
    range[0] = bounds[0];
    range[1] = bounds[1];
    // yaps_message(" start p(%lg)=%lg, range=[%lg,%lg]\n",x,y,range[0],range[1]);
    y += log(rng_unit(rng));
    // yaps_message("  y=%lg\n", y);
    for (tries=1; tries<TOOMANY; tries++) {
      x = range[0] + rng_unit(rng)*(range[1]-range[0]);
      if ( post(x,pars)>y ) {
	*xp = x;
	// yaps_message(" got p(%lg)=%lg after %d tries, range=[%lg,%lg]\n",
	//	    x, post(x,pars), tries, range[0], range[1]);
	break;
      }
      /*
       *  reduce bounds, assumes posterior is unimodal
       */
      if ( x<*xp ) 
	range[0] = x;
      else
	range[1] = x;
      // yaps_message("   p(%lg)=%lg, r=[%lg,%lg]\n",
      //             x,post(x,pars),range[0],range[1]);
    }
    if (tries>=TOOMANY ) {
      fprintf(stderr,
	      "SliceSimple: giving up after %d tries, range=[%lg,%lg]\n",
	      TOOMANY, range[0], range[1]);
      return 1;
    }
  }
  return 0;
}
