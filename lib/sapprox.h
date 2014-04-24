/*
 * Approximate S computation
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
 *  Exact/Approx calculations for S, for theory, see:
 *      http://arxiv.org/abs/1007.0296
 *
 *  Returns -HUGE_VAL if any trouble.
 *         e.g.,  m>4
 *          
 */

/*
 *    approximate calculation for M<=4, exact when a==0
 */
double S_approx(int n, int m, float a);

/*
 *    approx calculation for M<=4, derivative w.r.t. a
 */
double S_approx_da(int n, int m, float a) ;
