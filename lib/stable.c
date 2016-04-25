/*
 * Stirling Number table handling
 * Copyright (C) 2009-2014 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *   This is a rewrite of various earlier implementations with a different
 *   data layout and a few more options.
 *
 *   Note the code has a float and a double version.  The float version
 *   keeps two double vectors representating the current boundary of
 *   the Stirling number matrix.  This way it can be extended in either
 *   direction without sacrificing precision.  However, it makes
 *   extensions somewhat complicated, because one needs to keep
 *   track of the boundary during computation.
 *     
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "stable.h"
#include "yaps.h"
#ifdef S_USE_THREADS
#include <pthread.h>
#endif

#define MEMALLOCED
#ifdef MEMALLOCED
/********************************************************************
 *    this is a fudge to allow recording of memory allocation sizes;
 *    we store the allocation size along with the memory
 */
//   guess at basic memory requirements of malloc()
#define memsize(s) ((((s)+9)/8)*8)

 void *malloc_hook(stable_t *sp, size_t size) {
  void *ptr;
  size_t ms = memsize(size);
  size += sizeof (size_t);
  ptr = malloc(size);
  if ( !ptr ) return NULL;
  *(size_t *) ptr = size;
  if ( sp ) sp->memalloced += ms;
  return ((size_t *) ptr) + 1;
}

/*
 *    reallocs the memory at *ptr so that it is always
 *    available for use by anoter thread
 */
void realloc_hook(stable_t *sp, void **ptr, size_t size) {
  void *ptrtmp, *ptrsave;
  size_t ms = memsize(size);
  size_t oldsize = *(((size_t *)*ptr) - 1) - sizeof(size_t);
  assert(oldsize>0);
  if ( oldsize>=size )
 	return;
  size += sizeof(size_t);
  assert(oldsize<size);
  ptrtmp = malloc(size);
  if ( !ptrtmp ) {
    *ptr = NULL;
    return;
  }
  *(size_t *)ptrtmp = size;
  if ( sp ) sp->memalloced += ms - memsize(oldsize);
  ptrtmp = ((size_t *)ptrtmp) + 1;
  memcpy(ptrtmp, *ptr, oldsize);
  ptrsave = *ptr;
  *ptr = ptrtmp;
  free((void *)(((size_t *) ptrsave) - 1));
}

void free_hook (stable_t *sp, void *ptr) {
  if ( sp ) sp->memalloced -= memsize(* (((size_t *) ptr) - 1));
  ptr = (void *) (((size_t *) ptr) - 1);
  free(ptr);
}

#define myrealloc(p,x) realloc_hook(sp,(void**)&(p),x)
#define mymalloc(x) malloc_hook(sp,x)
#define myfree(x) free_hook(sp,x)
/********************************************************************/
#endif

static double logadd(double V, double lp) {
  if ( lp>V ) {
    // swap so V is bigger
    double t = lp;
    lp = V;
    V = t;
  }
  return V + log(1.0+exp(lp-V));
}

void S_tag(stable_t *S, char *tag) {
  S->tag = malloc(strlen(tag)+1);
  strcpy(S->tag,tag);
}

stable_t *S_make(unsigned initN, unsigned initM, unsigned maxN, unsigned maxM, 
 double a, uint32_t flags) {
  int N;
  stable_t *sp = NULL;
  sp = mymalloc(sizeof(stable_t));
  if ( !sp ) 
    return NULL;
  
  if ( maxM<10 )
    maxM = 10;
  if ( maxN<maxM )
    maxN = maxM;
  if ( initM<10 )
    initM = 10;
  if ( initN<initM )
    initN = initM;
  if ( initN>maxN )
    initN = maxM;
  if ( initN>maxN )
    initN = maxN;

  if ( (flags&S_STABLE)==0 && (flags&S_UVTABLE)==0 )
    return NULL;
  
  sp->tag = NULL;
  sp->memalloced = 0;
  sp->flags = flags;
  sp->maxN = maxN;
  sp->maxM = maxM;
  sp->usedN = initN;
  sp->usedM = initM;
  sp->usedN1 = initN;
  sp->startM = initM;
  sp->S = NULL;
  sp->SfrontN = sp->SfrontM = sp->S1 = NULL;
  sp->V = NULL;
  sp->VfrontN = sp->VfrontM = NULL;
  sp->Sf = sp->Vf = NULL;
#ifdef S_USE_THREADS
  if ( (flags&S_THREADS) ) {
    // yaps_message("Initialised mutex \n");
    pthread_mutex_init(&sp->mutex, NULL);
  }
#endif
  
  sp->S1 = mymalloc(sizeof(sp->S1[0])*(initN));
  if ( !sp->S1 ) {
    myfree(sp);
    return NULL;
  }
  if ( flags&S_STABLE ) {
    if ( flags&S_FLOAT ) {
      /*
       *  allocate frontier
       */
       sp->SfrontN = mymalloc(sizeof(sp->SfrontN[0])*(initM-1));
       if ( !sp->SfrontN ) {
         S_free(sp);  	return NULL;
       }
      /*
       *    sets diagonal entry of S since the loop writing
       *    SfrontN never does the diagnal itself
       */
       memset(sp->SfrontN,0,sizeof(sp->SfrontN[0])*(initM-1));
       sp->SfrontM = mymalloc(sizeof(sp->SfrontM[0])*(initN-initM+1));
       if ( !sp->SfrontM ) {
         S_free(sp);  	return NULL;
       }
      /*
       *  allocate sp->Sf[] as vector of vectors
       */
       sp->Sf = mymalloc(sizeof(sp->Sf[0])*(sp->usedN-2));
       if ( !sp->Sf ) {
         S_free(sp);  	return NULL;
       }
       memset(sp->Sf,0,sizeof(sp->Sf[0])*(sp->usedN-2));
      /*
       *   allocate sp->Sf[0][.] to sp->Sf[startM-3][.] in one block
       */
       sp->Sf[0] = mymalloc(sizeof(sp->Sf[0][0])*(sp->startM-1)*(sp->startM-2)/2);
       if ( !sp->Sf[0] ) {
         S_free(sp);  return NULL;
       }
       for (N=1; N<=sp->startM-3; N++) 
         sp->Sf[N] = sp->Sf[N-1] + N;
      /*
       *   allocate remaining sp->Sf[N][.] as vectors
       */
       assert(sp->startM-2+sp->Sf[sp->startM-3]-sp->Sf[0]==
        (sp->startM-1)*(sp->startM-2)/2);
       for (N=sp->startM-2; N<=sp->usedN-3; N++) {
         sp->Sf[N] = mymalloc(sizeof(sp->Sf[0][0])*(sp->usedM-1));
         if ( !sp->Sf[N] ) {
           S_free(sp);  return NULL;
         }
       }  
     } else { 
      sp->S = mymalloc(sizeof(sp->S[0])*sp->usedN);
      if ( !sp->S ) {
       S_free(sp);  return NULL;
     }
      /*
       *   allocate sp->S[0][.] to sp->S[startM-3][.] in one block
       */
       sp->S[0] = mymalloc(sizeof(sp->S[0][0])*(sp->startM-1)*(sp->startM-2)/2);
       if ( !sp->S[0] ) {
         S_free(sp);  return NULL;
       }
       for (N=1; N<=sp->startM-3; N++) 
         sp->S[N] = sp->S[N-1] + N;
      /*
       *   allocate remaining sp->S[N][.] as vectors for N>=usedM+1
       *   which store values M=2,...,usedM, so need (usedM-1) space
       */
       assert(sp->startM-2+sp->S[sp->startM-3]-sp->S[0]==
        (sp->startM-1)*(sp->startM-2)/2);
       for (N=sp->startM-2; N<=sp->usedN-3; N++) {
         sp->S[N] = mymalloc(sizeof(sp->S[0][0])*(sp->usedM-1));
         if ( !sp->S[N] ) {
           S_free(sp);  return NULL;
         }
       }  
     }
   }
   if ( flags&S_UVTABLE ) {
    if ( flags&S_FLOAT ) {
     /*
       *  allocate frontier
       */
       sp->VfrontN = mymalloc(sizeof(sp->VfrontN[0])*(initM-1));
       if ( !sp->VfrontN ) {
         S_free(sp);  	return NULL;
       }
      /*
       *    sets diagonal entry of V since the loop writing
       *    VfrontN never does the diagnal itself
       */
       memset(sp->VfrontN,0,sizeof(sp->VfrontN[0])*(initM-1));
       sp->VfrontM = mymalloc(sizeof(sp->VfrontM[0])*(initN-initM+1));
       if ( !sp->VfrontM ) {
         S_free(sp);  	return NULL;
       }
       
       sp->Vf = mymalloc(sizeof(sp->Vf[0])*sp->usedN);
       if ( !sp->Vf ) {
         S_free(sp);  return NULL;
       }
      /*
       *   allocate sp->Vf[0][.] to sp->Vf[startM-2][.] in one block
       */
       sp->Vf[0] = mymalloc(sizeof(sp->Vf[0][0])*(sp->startM-1)*(sp->startM)/2);
       if ( !sp->Vf[0] ) {
         S_free(sp);  return NULL;
       }
       for (N=1; N<=sp->startM-2; N++) 
         sp->Vf[N] = sp->Vf[N-1] + N;
      /*
       *   allocate remaining sp->Vf[N][.] as vectors for N>=usedM+1
       *   which store values M=2,...,usedM, so need (usedM-1) space
       */
       assert(sp->startM-1+sp->Vf[sp->startM-2]-sp->Vf[0]==
        (sp->startM-1)*(sp->startM)/2);
       for (N=sp->startM-1; N<=sp->usedN-2; N++) {
         sp->Vf[N] = mymalloc(sizeof(sp->Vf[0][0])*(sp->usedM-1));
         if ( !sp->Vf[N] ) {
           S_free(sp);  return NULL;
         }
       }  
     } else {
      sp->V = mymalloc(sizeof(sp->V[0])*sp->usedN);
      if ( !sp->V ) {
       S_free(sp);  return NULL;
     }
      /*
       *   allocate sp->V[0][.] to sp->V[startM-2][.] in one block
       */
       sp->V[0] = mymalloc(sizeof(sp->V[0][0])*(sp->startM-1)*(sp->startM)/2);
       if ( !sp->V[0] ) {
         S_free(sp);  return NULL;
       }
       for (N=1; N<=sp->startM-2; N++) 
         sp->V[N] = sp->V[N-1] + N;
      /*
       *   allocate remaining sp->V[N][.] as vectors for N>=usedM+1
       *   which store values M=2,...,usedM, so need (usedM-1) space
       */
       assert(sp->startM-1+sp->V[sp->startM-2]-sp->V[0]==
        (sp->startM-1)*(sp->startM)/2);
       for (N=sp->startM-1; N<=sp->usedN-2; N++) {
         sp->V[N] = mymalloc(sizeof(sp->V[0][0])*(sp->usedM-1));
         if ( !sp->V[N] ) {
           S_free(sp);  return NULL;
         }
       }  
     }
   }

  /*
   *  this is where we actually build the Stirling numbers
   */
   S_remake(sp,a);
   return sp;
 }


/*
 *    assumes s->usedN/M already set to new values and memory filled
 *    startN/M = 0  --->  refill everything
 *    startN/M > 0  --->  memory extended so refill from here,
 *                        i.e.,  these were *last* values set, start +1
 */
 static int S_remake_part(stable_t *sp, double a, 
			  unsigned startN, unsigned startM,
			  unsigned usedN, unsigned usedM, unsigned usedN1) {
  int N, M;
  
  if ( startN==0 )
    startM = 0;
  sp->a = a;
  sp->lga = lgamma(1.0-a);

  //  yaps_message("S_remake_part(a=%lf, start: N=%u, M=%u, used: N=%u, M=%u, usedN1=%u)\n", a, startN, startM, usedN, usedM, usedN1);

  /*
   *  need to reset at sp->S1[] least to usedN1;
   *  up to usedN used by sp->S[][],
   *  and data needs overwriting up to usedN1
   */
   if ( startN==0 ) {
    sp->S1[0] = 0;
    N = 2;
  } else {
    N = startN+1;
    // assert(sp->S1[startN-1]>0); ???
    if ( sp->S1[startN-1]==0 )
      sp->S1[startN-1] = lgamma(startN - sp->a) - sp->lga;
  }
  for ( ; N<=usedN; N++)
    sp->S1[N-1] = sp->S1[N-2] + log(N-1-a);
  
  if ( startN==0 )
    //   a has changed, so reset others
    for ( ; N<=usedN1; N++)
      sp->S1[N-1] = 0;
  
  if ( sp->flags&S_STABLE ) {
    if ( (sp->flags&S_FLOAT)==0 ) {
      if ( startM>0 && startM<usedM ) {
	/*
	 *   extend for M upto startN
	 */
	for (N=startM+1; N<=startN; N++) {
	  for (M=startM+1; M<N && M<=usedM; M++) {
	    sp->S[N-3][M-2] = 
	      logadd(log(N-M*a-1.0)+((M<N-1)?sp->S[N-4][M-2]:0), 
		     sp->S[N-4][M-3]);
	    assert(ISFINITE(sp->S[N-3][M-2]));
	  }
	}
      }
      /*
       *   now fill from 0 after startN ... just like usual
       */
      if ( startN==0 ) {
	sp->S[0][0] = logadd(sp->S1[1],log(2-2*a));
	N = 4;
      } else {
	N = startN+1;
	assert(sp->S[N-4][0]>0);
      }
      for (; N<=usedN; N++) {
	sp->S[N-3][0] = logadd(log(N-2*a-1.0)+sp->S[N-4][0], 
			       sp->S1[N-2]);
	for (M=3; M<=usedM && M<N; M++) {
	  sp->S[N-3][M-2] = 
	    logadd(log(N-M*a-1.0)+((M<N-1)?sp->S[N-4][M-2]:0), sp->S[N-4][M-3]);
	  assert(ISFINITE(sp->S[N-3][M-2]));
	}
       }
    } else {
      /*
       *   computation done in double by storing in sp->SfrontN+M[]
       */
      if ( startM>0 && startM<usedM ) {
	/*
	 *   extend for M upto startN
	 */
	for (N=startM+2; N<=startN; N++) {
	  double lastS;
	  if (startM+1<N-1)
	    lastS = sp->Sf[N-4][startM-1];
	  else
	    lastS = 0 ;
	  sp->Sf[N-3][startM-1] = sp->SfrontN[startM-1]
	    = logadd(log(N-(startM+1)*a-1.0)+lastS,
		     sp->SfrontM[N-startM-2]);
	  for (M=startM+2; M<N && M<=usedM; M++) {
	    double saveS = sp->SfrontN[M-2];
	    if ( M==N-1 ) saveS = 0;
	    sp->SfrontN[M-2] = logadd(log(N-M*a-1.0)+saveS, lastS);
	    sp->Sf[N-3][M-2] = sp->SfrontN[M-2];
	    assert(ISFINITE(sp->Sf[N-3][M-2]));
	    lastS = saveS;
	  }
	  //   save the SfrontM value
	  if ( N>usedM )
	    sp->SfrontM[N-usedM-1] = sp->SfrontN[usedM-2];
	}
      }
      if ( startN==0 ) {
	sp->Sf[0][0] = sp->SfrontN[0] = logadd(sp->S1[1],log(2-2*a));
	N = 4;
      } else {
	N = startN+1;
	assert(sp->Sf[N-4][0]>0);
      }
      for ( ; N<=usedN; N++) {
	double lastS;
	lastS = sp->SfrontN[0];
	sp->Sf[N-3][0] = sp->SfrontN[0] =
	  logadd(log(N-2*a-1.0)+lastS, sp->S1[N-2]);
	for (M=3; M<=usedM && M<N; M++) {
	  double saveS = sp->SfrontN[M-2];
	  if (M==N-1) saveS = 0;
	  sp->SfrontN[M-2] =
	    logadd(log(N-M*a-1.0)+saveS, lastS);
	  sp->Sf[N-3][M-2] = sp->SfrontN[M-2];
#ifndef NDEBUG
	  if ( !ISFINITE(sp->Sf[N-3][M-2]) ) 
	    yaps_quit("Building '%s' N to %d, sp->Sf[%d][%d] not finite, from %lf,%lf\n",
		      sp->tag, usedN, N-3, M-2,  saveS, lastS);
#endif
	  assert(ISFINITE(sp->Sf[N-3][M-2]));
	  lastS = saveS;
	}
	//   save the SfrontM value
	if ( N>usedM )
	  sp->SfrontM[N-usedM-1] = sp->SfrontN[usedM-2];
      }
    }
  }    
  if ( sp->flags&S_UVTABLE ) {
    if ( (sp->flags&S_FLOAT)==0 ) {
      if ( startM>0 && startM<usedM ) {
	/*
	 *   extend for M upto startN
	 */
	for (N=startM+1; N<=startN; N++) {
	  for (M=startM+1; M<=N && M<=usedM; M++) {
	    sp->V[N-2][M-2] = 
	      (1.0+((M<N)?((N-1-M*a)*sp->V[N-3][M-2]):0))
	      / (1.0/sp->V[N-3][M-3]+(N-1-(M-1)*a));
	  }
	}
      }
      /*
       *   now fill from 0 after startN ... just like usual
       */
      if ( startN==0 ) {
	sp->V[0][0] = 1.0/(1.0-a);
	N = 3;
      } else {
	N = startN+1;
	assert(sp->V[N-3][0]>0);
      }
      for (; N<=usedN; N++) {
	sp->V[N-2][0] = (1.0+(N-1-2*a)*sp->V[N-3][0])/(N-1-a);
	for (M=3; M<=usedM && M<=N; M++) {
	  sp->V[N-2][M-2] = 
	    (1.0+((M<N)?((N-1-M*a)*sp->V[N-3][M-2]):0))
	    / (1.0/sp->V[N-3][M-3]+(N-1-(M-1)*a));
	}
      }
    } else {
      if ( startM>0 && startM<usedM ) {
	/*
	 *   extend for M upto startN
	 */
	for (N=startM+1; N<=startN; N++) {
	  double lastS;
	  if ( startM+1<N) 
	    lastS = sp->VfrontN[startM-1];
	  else 
	    lastS = 0;  
	  sp->Vf[N-2][startM-1] = sp->VfrontN[startM-1] = 
	    (1.0+(N-1-(startM+1)*a)*lastS)
	    / (1.0/sp->VfrontM[N-1-startM]+(N-1-(startM)*a));
	  for (M=startM+2; M<=N && M<=usedM; M++) {
	    double saveS = sp->VfrontN[M-2];
	    assert(lastS!=0);
	    sp->Vf[N-2][M-2] = sp->VfrontN[M-2] = 
	      (1.0+((M<N)?((N-1-M*a)*saveS):0))
	      / (1.0/lastS+(N-1-(M-1)*a));
	    lastS = saveS;
	  }
	  //   save the VfrontM value
	  if ( N>=usedM )
	    sp->VfrontM[N-usedM] = sp->VfrontN[usedM-2];
	}
      }
      /*
       *   now fill from 0 after startN ... just like usual
       */
      if ( startN==0 ) {
	sp->Vf[0][0] = sp->VfrontN[0] = 1.0/(1.0-a);
	N = 3;
      } else {
	N = startN+1;
	assert(sp->Vf[N-3][0]>0);
      }
      for (; N<=usedN; N++) {
	double lastS;
	lastS = sp->VfrontN[0];
	sp->Vf[N-2][0] = sp->VfrontN[0] =
	  (1.0+(N-1-2*a)*lastS)/(N-1-a);
	for (M=3; M<=usedM && M<=N; M++) {
	  double saveS = sp->VfrontN[M-2];
	  assert(lastS!=0);
	  sp->Vf[N-2][M-2] = sp->VfrontN[M-2] = 
	    (1.0+((M<N)?((N-1-M*a)*saveS):0))
	    / (1.0/lastS+(N-1-(M-1)*a));
	  lastS = saveS;
	}
	//   save the VfrontM value
	if ( N>=usedM )
	  sp->VfrontM[N-usedM] = sp->VfrontN[usedM-2];
      }
    }
  }
  /*
   *  change bounds at end only after data filled;
   *  in case other threads running
   */
  sp->usedN = usedN;
  sp->usedN1 = usedN1;
  sp->usedM = usedM;
  return 0;
 }

 int S_remake(stable_t *sp, double a) {
   int ret = S_remake_part(sp, a, 0, 0, sp->usedN, sp->usedM, sp->usedN1);
  if ( !ret && (sp->flags&S_VERBOSE) ) 
    S_report(sp, stderr);
  return ret;
}

/*
 *    assume tables filled to usedN,usedM;
 *    request (N,M) table value;
 *    fiddle maxM, maxN to make it non-trivial
 *    create extra space first;
 *    then remake table values;
 *    return non-zero on error
 */
 static int S_extend(stable_t *sp, int N, int M) {
  int n;
  int result = 0;
  //  int Nin = N, Min = M;
  N++; M++;
  /*
   *  N shouldn't be too big or small
   */
#ifdef S_USE_THREADS
   if ( (sp->flags&S_THREADS) ) 
    pthread_mutex_lock(&sp->mutex);
#endif
  if ( N<sp->usedN && M<sp->usedM  ) 
    /*
     *  someone did the change before we got here!
     */
   goto SE_result;
   
   if ( N<sp->usedN )
    N = sp->usedN;
  if ( N>sp->maxN )
    N = sp->maxN;

  /*
   *   N increase should not be trivial
   */
   if ( N>sp->usedN ) {
    if ( N<sp->usedN*1.1 )
      N = sp->usedN*1.1;
    if ( N<sp->usedN+50 )
      N = sp->usedN+50;
    /*
     *  reset if made too big
     */
     if ( N>sp->maxN ) {
      N = sp->maxN;
    }
  }

  /*
   *  M shouldn't be too big or small
   */
   if ( M<sp->usedM )
    M = sp->usedM;
  if ( N<M )
    M = N;
  if ( M>sp->maxM )
    M = sp->maxM;

  /*
   *   M increase should not be trivial
   */
   if ( M>sp->usedM ) {
    if ( M<sp->usedM*1.1 )
      M = sp->usedM*1.1;
    if ( M<sp->usedM+50 )
      M = sp->usedM+50;
    /*
     *  reset if made too big
     */
     if ( M>sp->maxM ) {
      M = sp->maxM;
    }
    if ( M>sp->usedN ) {
      M = sp->usedN;
    }
  }
  /*
   *   N and M values now set
   */

   if ( M>sp->usedM ) {
    /*
     *    extend size of existing vectors;
     *    note S/Sf/V/Vf are triangular up to n=usedM, so ignore that
     */
     for (n=sp->usedM+1; n<=sp->usedN; n++) {
      if ( sp->flags&S_STABLE ) {
       if ( (sp->flags&S_FLOAT)==0 ) {
         myrealloc(sp->S[n-3],sizeof(sp->S[0][0])*(M-1));
         if ( !sp->S[n-3] ) {
           S_free(sp);
           result = 1;
           goto SE_result;
         }
       } else {
         myrealloc(sp->Sf[n-3],sizeof(sp->Sf[0][0])*(M-1));
         if ( !sp->Sf[n-3] ) {
           S_free(sp);
           result = 1;
           goto SE_result;
         }
       }
     }
     if ( sp->flags&S_UVTABLE ) {
       if ( (sp->flags&S_FLOAT)==0 ) {
         myrealloc(sp->V[n-2],sizeof(sp->V[0][0])*(M-1));
         if ( !sp->V[n-2] ) {
           S_free(sp);
           result = 1;
           goto SE_result;
         }
       } else {
         myrealloc(sp->Vf[n-2],sizeof(sp->Vf[0][0])*(M-1));
         if ( !sp->Vf[n-2] ) {
           S_free(sp);
           result = 1;
           goto SE_result;
         }
       }
     }
     }
     if ( sp->flags&S_STABLE && sp->flags&S_FLOAT ) {
       myrealloc(sp->SfrontN,sizeof(sp->SfrontN[0])*(M-1));
       if ( !sp->SfrontN ) {
	 S_free(sp);
	 result = 1;
	 goto SE_result;
       }
     }
     if ( sp->flags&S_UVTABLE && sp->flags&S_FLOAT ) {
       myrealloc(sp->VfrontN,sizeof(sp->VfrontN[0])*(M-1));
       if ( !sp->VfrontN ) {
	 S_free(sp);
	 result = 1;
	 goto SE_result;
       }
     }
   }
   /*
    *    now create new vectors in S/Sf/V/Vf
    */
   if ( N>sp->usedN ) {
     /*
      *   extend size of S/Sf/SfrontN
      */
     if ( sp->flags&S_STABLE ) {
       myrealloc(sp->S1,sizeof(sp->S1[0])*N);
       if ( !sp->S1 ) {
	 S_free(sp);
	 result = 1;
	 goto SE_result;
       }  
       if ( (sp->flags&S_FLOAT)==0 ) {
	 myrealloc(sp->S,sizeof(sp->S[0])*N);
	 if ( !sp->S ) {
	   S_free(sp);
	   result = 1;
	   goto SE_result;
	 }
       } else {
	 myrealloc(sp->Sf,sizeof(sp->Sf[0])*N);
	 if ( !sp->Sf ) {
	   S_free(sp);
	   result = 1;
	   goto SE_result;
	 }
	 myrealloc(sp->SfrontM,sizeof(sp->SfrontM[0])*(N-M));
	 if ( !sp->SfrontM ) {
	   S_free(sp);
	   result = 1;
	   goto SE_result;
	 }
       }
     }
   }
   /*
    *    extend size of V/Vf/VfrontN+M
    */
   if ( sp->flags&S_UVTABLE ) {
     if ( (sp->flags&S_FLOAT)==0 ) {
       myrealloc(sp->V,sizeof(sp->V[0])*N);
       if ( !sp->V ) {
	 S_free(sp);
	 result = 1;
	 goto SE_result;
       }
     } else {
       myrealloc(sp->Vf,sizeof(sp->Vf[0])*N);
       if ( !sp->Vf ) {
	 S_free(sp);
	 result = 1;
	 goto SE_result;
       }
       if ( N-M > sp->usedN-sp->usedM ) {
	 /*
	  *   stores *last* values, so cannot shrink or loose a few
	  */
	 myrealloc(sp->VfrontM,sizeof(sp->VfrontM[0])*(N-M+1));
	 if ( !sp->VfrontM ) {
	   S_free(sp);
	   result = 1;
	   goto SE_result;
	 }
       }
     }
   }
   /*
    *    now create new vectors
    */
   for (n=sp->usedN+1; n<=N; n++) {
     if ( sp->flags&S_STABLE ) {
       if ( (sp->flags&S_FLOAT)==0 ) {
	 sp->S[n-3] = mymalloc(sizeof(sp->S[0][0])*(M-1));
       if ( !sp->S[n-3] ) {
         S_free(sp);
         result = 1;
         goto SE_result;
       }
       } else {
	 sp->Sf[n-3] = mymalloc(sizeof(sp->Sf[0][0])*(M-1));
	 if ( !sp->Sf[n-3] ) {
	   S_free(sp);
	   result = 1;
	   goto SE_result;
	 }
       }
     }
     if ( sp->flags&S_UVTABLE ) {
       if ( (sp->flags&S_FLOAT)==0 ) {
	 sp->V[n-2] = mymalloc(sizeof(sp->V[0][0])*(M-1));
	 if ( !sp->V[n-2] ) {
	   S_free(sp);
	   result = 1;
	   goto SE_result;
	 }
       } else {
	 sp->Vf[n-2] = mymalloc(sizeof(sp->Vf[0][0])*(M-1));
	 if ( !sp->Vf[n-2] ) {
	   S_free(sp);
	   result = 1;
	   goto SE_result;
	 }
       }
     }
   }
   {
     int oldN, oldM;
     oldN = sp->usedN;
     oldM = sp->usedM;
     result = S_remake_part(sp,sp->a, oldN, oldM, N, M, sp->usedN1);
   }
 SE_result:
#ifdef S_USE_THREADS
   if ( (sp->flags&S_THREADS) ) {
     // yaps_message("Extended %s under lock:  in=%d,%d  set=%d,%d  used=%d,%d\n",
		   // sp->tag, Nin, Min, N, M, sp->usedN, sp->usedM);
     pthread_mutex_unlock(&sp->mutex);
   }
#endif
   return result;
 }

/*
 *   using cache if possible,
 *   has its own maximum memory allowed to be
 *   greater than usedN
 */
double S_S1(stable_t *sp, unsigned n) {
  // int nin = n;
  if ( n==0 )
    return -HUGE_VAL;
  if ( !sp->S1 )
    return -HUGE_VAL;
  if ( n>sp->usedN ) {
#ifdef S_USE_THREADS
    if ( (sp->flags&S_THREADS) ) 
      pthread_mutex_lock(&sp->mutex);
#endif
    if ( n>sp->usedN ) {
      /*
       *   value may not be set here
       */
      if ( n>sp->usedN1 ) {
	/*
	 *   possibly extend memory and initialise
	 */
	int i;
	if ( n>sp->maxN && !(sp->flags & S_ASYMPT) ) {
	  return -HUGE_VAL;
	}
	if ( n<sp->usedN1*1.1 )
	  n = sp->usedN1*1.1;
	if ( n<sp->usedN1+50 )
	  n = sp->usedN1 + 50;
	if ( n>sp->maxN )
	  n = sp->maxN;
	myrealloc(sp->S1, sizeof(sp->S1[0])*n);
	if ( !sp->S1 ) {
	  return -HUGE_VAL;
	}
	for (i=sp->usedN1; i<n; i++) 
	  sp->S1[i] = 0;
	sp->usedN1 = n; 
      }
      if ( sp->S1[n-1]==0 ) {
	if ( sp->S1[n-2]==0 )
	  sp->S1[n-1] = lgamma(n-sp->a) - sp->lga;
	else
	  sp->S1[n-1] = sp->S1[n-2] + log(n-1-sp->a);
      }
    }
    // yaps_message("S_S1(%d): usedN1=%d, usedN=%d\n", nin, sp->usedN1, sp->usedN);
#ifdef S_USE_THREADS
    if ( (sp->flags&S_THREADS) ) 
      pthread_mutex_unlock(&sp->mutex);
#endif
  }
  return sp->S1[n-1];
}

double S_U(stable_t *sp, unsigned n, unsigned m) {
  if ( m==1 )
    return n - sp->a;
  if ( m<=1 ) 
    yaps_quit("Bad constraints in S_U(%s,%u,%u)\n",
      sp->tag, n, m);
  assert(m>1);
  return n - m*sp->a + 1/S_V(sp,n,m);
}

double S_UV(stable_t *sp, unsigned n, unsigned m) {
  double SV;
  if ( m==1 )
    return -HUGE_VAL;
  /*  identity because S^n_n==1 */
  if ( m==n+1 )
   return 1;
  /*  identity because U^n_n = n(n+1)(1-a)/2  */
  if ( m==n )
    return (n+1.0)/(n-1.0);
  SV = S_V(sp,n,m); 
  return (n - m*sp->a)*SV + 1.0;
}


double S_V(stable_t *sp, unsigned n, unsigned m) {
  if ( (sp->flags & S_UVTABLE)==0 )
    return 0;
  if ( m>=sp->usedM-1 || n>=sp->usedN-1 ) {
    if ( n>sp->maxN || m>sp->maxM  ) {
      if ( n>sp->maxN && (sp->flags & S_ASYMPT) ) {
	if ( sp->a>0 )
	  return (1.0-pow(n,-sp->a))/sp->a/(m-1);
	else {
	  double ln = log(n);
	  return  ln/(m-1) * exp(gamma(1+(m-2)/ln)-gamma(1+(m-1)/ln));
	}
      }	  
      if ( (sp->flags & S_QUITONBOUND) ) {
	assert(n>sp->maxN || m>sp->maxM);
	if ( sp->tag )
	  yaps_quit("S_V(%u,%u,%lf) tagged '%s' hit bounds (%u,%u)\n",n,m,sp->a,
		    sp->tag, sp->maxN, sp->maxM);
	else 
	  yaps_quit("S_V(%u,%u,%lf) hit bounds\n",n,m,sp->a);
      } else
	return 0;
    }
    // yaps_message("S_V(%s,%d,%d) calling extend\n", sp->tag, n, m);
    if ( S_extend(sp,n+1,m+1) ) 
      yaps_quit("S_extend() out of memory\n");
  }
  assert(m>=2);
  if ( n<m ) return 0;
  if ( sp->flags & S_FLOAT ) {
    assert(sp->Vf);
    if ( sp->Vf[n-2]==NULL )
      yaps_quit("S_V(%s,%u,%u) Vf memory unavailable\n", sp->tag, n, m);
    assert(sp->Vf[n-2]);
    return sp->Vf[n-2][m-2];
  }
  assert(sp->V);
  assert(sp->V[n-2]);
  return sp->V[n-2][m-2];
}

double S_S(stable_t *sp, unsigned N, unsigned T) {
  if ( (sp->flags & S_STABLE)==0 )
    return -HUGE_VAL;
  if ( N==T )
    return 0;
  if ( T==1 )
    return S_S1(sp, N);
  if ( N<T || T==0 )
    return -HUGE_VAL;
  if ( T>sp->usedM || N>sp->usedN ) {
    if ( N>sp->maxN || T>sp->maxM  ) {
      if ( N>sp->maxN && (sp->flags & S_ASYMPT) )
	return S_asympt(sp, N, T);
      if ( (sp->flags & S_QUITONBOUND) )
       if ( sp->tag )
         yaps_quit("S_S(%u,%u,%lf) tagged '%s' hit bounds\n",N,T,sp->a,
          sp->tag);
       else 
         yaps_quit("S_S(%u,%u,%lf) hit bounds\n",N,T,sp->a);
       else
         return -HUGE_VAL;
     }
     if ( S_extend(sp,N+1,T+1) )
      yaps_quit("S_extend() out of memory\n");
  }
  if ( sp->flags&S_FLOAT ) {
    assert(sp->Sf);
    assert(sp->Sf[N-3]);
    return sp->Sf[N-3][T-2];
  } 
  assert(sp->S);
  assert(sp->S[N-3]);
  return sp->S[N-3][T-2];
}

/*
 *    only frees allocated memory, can be called
 *    during a failed make too
 */
 void S_free(stable_t *sp) {
  if ( !sp )
    return;
  if ( sp->tag ) free(sp->tag);
  if ( sp->S1 )  myfree(sp->S1);
  if ( sp->SfrontN )  myfree(sp->SfrontN);
  if ( sp->SfrontM )  myfree(sp->SfrontM);
  if ( sp->VfrontN )  myfree(sp->VfrontN);
  if ( sp->VfrontM )  myfree(sp->VfrontM);
  if ( sp->S ) {
    int n;
    for (n=sp->startM-2; n<=sp->usedN-3; n++)
      if ( sp->S[n] ) myfree(sp->S[n]);
    myfree(sp->S[0]);
    myfree(sp->S);
  }
  if ( sp->Sf ) {
    int n;
    for (n=sp->startM-2; n<=sp->usedN-3; n++)
      if ( sp->Sf[n] ) myfree(sp->Sf[n]);
    myfree(sp->Sf[0]);
    myfree(sp->Sf);
  }
  if ( sp->V ) {
    int n;
    for (n=sp->startM-1; n<=sp->usedN-2; n++)
      if ( sp->V[n] ) myfree(sp->V[n]);
    myfree(sp->V[0]);
    myfree(sp->V);
  }
  if ( sp->Vf ) {
    int n;
    for (n=sp->startM-1; n<=sp->usedN-2; n++)
      if ( sp->Vf[n] ) myfree(sp->Vf[n]);
    myfree(sp->Vf[0]);
    myfree(sp->Vf);
  }
#ifdef S_USE_THREADS
  if ( (sp->flags&S_THREADS) ) {
    pthread_mutex_destroy(&sp->mutex);
  }
#endif
  myfree(sp);
}

void S_report(stable_t *sp, FILE *fp) {
  if ( fp ) {
    if ( sp->tag ) 
      fprintf(fp, "S-table '%s': ", sp->tag);
    else
      fprintf(fp, "S-table: ");
    fprintf(fp, "a=%lf, N=%u/%u, M=%u/%u, %s%s %s",
     sp->a, sp->usedN, sp->maxN, sp->usedM, sp->maxM, 
     (sp->flags&S_STABLE)?"+S":"", 
     (sp->flags&S_UVTABLE)?"+U/V":"", 
     (sp->flags&S_FLOAT)?"float":"double");
#ifdef MEMALLOCED
    fprintf(fp, " mem=%uk\n", sp->memalloced/1024);
#endif
    fprintf(fp, "\n");
  } else {
    if ( sp->tag ) 
      yaps_message("S-table '%s': ", sp->tag);
    else
      yaps_message("S-table: ");
    yaps_message("a=%lf, N=%u/%u, M=%u/%u, %s%s %s",
     sp->a, sp->usedN, sp->maxN, sp->usedM, sp->maxM, 
     (sp->flags&S_STABLE)?"+S":"", 
     (sp->flags&S_UVTABLE)?"+U/V":"", 
     (sp->flags&S_FLOAT)?"float":"double");
#ifdef MEMALLOCED
    yaps_message(" mem=%uk", sp->memalloced/1024);
#endif
    yaps_message("\n");
  }
} 

double S_asympt(stable_t *sp, unsigned n, unsigned m) {
  if ( sp->a==0 ) {
    /*
     *   https://www.researchgate.net/publication/2415504_Asymptotic_Expansions_for_the_Stirling_Numbers_of_the_First_Kind
     *   Asymptotic Expansions for the Stirling Numbers of the First Kind
     *   Hsien-Kuei Hwang, 2001
     */
    double ln = log(n);
    return gamma(n) + (m-1)*log(ln) - gamma(m) - gamma(1+(m-1)/ln);
  } else {
    double prod = 0;
    double la1 = lgamma(1.0-sp->a);
    double aln = sp->a*log((double)n);
    double np = pow(n,-sp->a);
    /*
     *    this works pretty well, but want a small factor
     *    of pow(n, ...), say  pow(0.2*n, ...)
     *    make that forma table
     */
    prod += lgamma((double)n) - la1 - lgamma((double)m)
      - (m-1.0)*log(sp->a) - aln;
    if ( np<1e-5 ) {
      prod -= (m-1)*np*(1+np*(0.5+np/3.0));
    } else
      prod += (m-1)*log(1.0-np);
    return prod;
  }
}
