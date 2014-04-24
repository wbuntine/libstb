/*
 * digamma routines
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
 */

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "digamma.h"

#ifndef LS_NOPOLYGAMMA
/*
 *    Minka's algorithm
 *    gets overflow on values > 80;
 *    large negative values return near zero;
 */
double digammaInv(double x) {
  double guess;
  int i;
  if ( x < -2.22 )
    guess = -1 / (x - digamma(1.0));
  else
    guess = exp(x) + 0.5;
  for (i=0; i<5; i++ ) {
    guess -= (digamma(guess) - x)/ trigamma(guess);
  }
  return guess;
}
#endif
