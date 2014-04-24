/*
 * caching difference of lgamma and digamma
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
 *     
 */

#ifndef __LGAMMA_H
#define __LGAMMA_H

/*
 *   caches are self allocating
 */
#define GCACHE 100
struct gcache_s {
  double par;
  double lgpar;
  double cache[GCACHE];
} ;

void gcache_init(struct gcache_s *lpg, double p);
double gcache_value(struct gcache_s *lpg, int j);
void pcache_init(struct gcache_s *lpg, double p);
double pcache_value(struct gcache_s *lpg, int j);
void qcache_init(struct gcache_s *lpg, double p);
double qcache_value(struct gcache_s *lpg, int j);

double gammadiff(int N, double alpha, double lga);
double psidiff(int N, double alpha, double pa);

#endif

