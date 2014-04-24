/*
 * List values for different Stirling number functions and ratios
 * Copyright (C) 2011-2012 Wray Buntine and Lan Du
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *         Lan Du (lan.du@nicta.com.au)
 *     
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yaps.h"
#include "stable.h"

/*
 *   parameters and hyper parameters
 */
float  apar = 0.0, bpar = 1.0;

int N=200;
int T=30;

void usage() {
  fprintf(stderr,"Commandline:  [-a val][-N val][-T val]\n");
  fprintf(stderr,"    print table of values for S_V(), S_U() and S_S()\n");
  fprintf(stderr,
          "Options are:\n"
	  "   -a val    # discount par (default=%f)\n"
	  "   -h           # print help message\n"
	  "   -N samples   # number of samples (default=%d)\n"
	  "   -T max       # maximum t for building tables (default=%d)\n",
	  apar, N, T);
}

int main(int argc, char* argv[])
{
  int c;
  int n, t;
  stable_t *S;

  /*
   *  default values for args
   */
  while ( (c=getopt(argc, argv,"ha:N:T:"))>=0 ) {
    switch ( c ) {
    case 'a':
      if ( !optarg || sscanf(optarg,"%f",&apar)!=1 )
	yaps_quit("Need a valid 'a' argument\n");
      break;
    case 'N':
      if ( !optarg || sscanf(optarg,"%d",&N)!=1 )
	yaps_quit("Need a valid 'N' argument\n");
      break;
    case 'T':
      if ( !optarg || sscanf(optarg,"%d",&T)!=1 )
	yaps_quit("Need a valid 'T' argument\n");
      break;
    case 'h':
      usage();
      exit(0);
    }
  }
  if ( T>N ) T = N;
      
  if ( !(S=S_make(N+1,T,2*N,2*T,apar,S_UVTABLE|S_STABLE|S_VERBOSE)) )
    yaps_quit("S_make failed\n");
  
  S_report(S, stdout);

#if 1
  /*
   *   list various values
   */
  printf("S(%d,%d) = %10.6lg, V=n/a U=%lg\n", N, 1, S_S(S,N,1), S_U(S,N,1));
  for (t=2; t<=T; t++) 
    printf("S(%d,%d) = %10.6lg, V=%lg U=%lg\n", N, t, 
	   S_S(S,N,t), S_V(S,N,t), S_U(S,N,t));

  printf("\nS(%d,%d) = %10.6lg, V=n/a U=%lg\n", 
	 N+10, 1, S_S(S,N+10,1), S_U(S,N+10,1));
  for (t=2; t<=2*T+10; t++) 
    printf("S(%d,%d) = %10.6lg, V=%lg U=%lg\n", N+10, t, 
	   S_S(S,N+10,t), S_V(S,N+10,t), S_U(S,N+10,t));
  for (n=T+1; n<=2*N; n++) 
    printf("S(%d,%d..) = %10.6lg %10.6lg\n", n, T, S_S(S,n,T), S_S(S,n,T+1));
  
  for (n=2; n<8; n++ ) {
    printf("S(%d,%d) = %10.6lg ", n, 1, S_S(S,n,1));
    for (t=2; t<=n; t++)
      printf(" %10.6lg", S_S(S,n,t));
    printf("\n");
  }
  for (n=2; n<8; n++ ) {
    printf("V(%d,%d) = %10.6lg", n, 2, S_V(S,n,2));
    for (t=3; t<=n; t++)
      printf(" %10.6lg", S_V(S,n,t));
    printf("\n");
  }
  for (n=2; n<8; n++ ) {
    printf("U(%d,%d) = %10.6lg", n, 2, S_U(S,n,2));
    for (t=3; t<=n; t++)
      printf(" %10.6lg", S_U(S,n,t));
    printf("\n");
  }
#else
  {
    /*
     *     Sample a partition of size T of N by Chinese Rest. distribution
     *     start by sampling the last entry from (1,2,...,N-T+1);
     *     see p(m | CRD, apr, N, T) on page 4 of "doc/alpha.pdf"
     */
    double *prob = malloc(sizeof(*prob)*N);
    double probtot = 1.0;
    prob[1] = 1.0;
    for (t=2; t<=N-T+1; t++) 
      probtot += prob[t] = (N-t+1)*(t-apar) / S_U(S,N-t,T-1)/(t-1) * prob[t-1]; 
    for (t=1; t<=N-T+1; t++) 
      printf("p(m=%d) = %lg\n", t, prob[t]/probtot );
  }
#endif

  S_report(S, stdout);
  S_free(S);
  return 0; 
}
