/*
 *   This is a simple test.
 *   (C) Wray Buntine 2011
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 *   It does side-by-side float and double calculations of
 *   the log computation and the ratio method for Stirling numbers.
 *   The ratio recursion is far more accurate.
 *
 *   Compile with:
 *        cc -o pt precision_test.c -lm -lrt
 *   then run
 *        ./pt
 *
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p)
{
  return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
           ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}

double logadd(double V, double lp) {
  if ( lp>V ) {
    // swap so V is bigger
    double t = lp;
    lp = V;
    V = t;
  }
  return V + log(1.0+exp(lp-V));
}

float flogadd(float V, float lp) {
  if ( lp>V ) {
    // swap so V is bigger
    float t = lp;
    lp = V;
    V = t;
  }
  return V + log(1.0+exp(lp-V));
}

#ifdef try_TS
double **TStbl=NULL;
static int TS_Nmax, TS_Mmax;
static double TS_a, TSlga;
static double *TS1 = NULL;

/*
 *   return log S^n_{1,a} = lgamma(n-a) - lgamma(1-a)
 *   using cache if possible
 */
double TS_S1(int n) {
  assert(n>0);
  assert(n<=TS_Nmax);
  if ( n==1 )
    return 0;
  if ( !isfinite(TS1[n-1]) ) {
    if ( !isfinite(TS1[n-2]) )
      TS1[n-1] = lgamma(n-TS_a) - TSlga;
    else
      TS1[n-1] = TS1[n-2] + log(n-1-TS_a);
  }
  return TS1[n-1];
}

double TS_safe(int n, int m) {
  if (m==0) {
    if ( n==0 ) return 1;
    return 0;
  } else if ( m==1 )
    return TS_S1(n);
  else
    return TS_S1(n)  + log(TStbl[m-2][n-1]);
}

double TS_make(int Nmax, int Mmax, double a) {
  int n,m;
  TS_a = a;
  TSlga = lgamma(1-a);
  TS1 = calloc(Nmax,sizeof(*TS1));
  TS1[0] = 0;
  TS1[1] = log(1.0-a);
  TS1[2] = TS1[1] + log(2.0-a);
  for (n=3; n<Nmax; n++) TS1[n] = -INFINITY;

  TStbl = malloc((Mmax-1)*sizeof(*TStbl));
  TStbl[0] = malloc((Mmax-1)*(Nmax-1)*sizeof(**TStbl));
  for (m=2; m<Mmax; m++)
    TStbl[m-1] = TStbl[m-2] + Nmax-1;

  TS_Nmax = Nmax;
  TS_Mmax = Mmax;
  TStbl[0][1] = 1.0/(1.0-a);
  for (n=2; n<Nmax; n++) {
    TStbl[0][n] =  (1.0 + (n-2.0*a)*TStbl[0][n-1]) / (n-a);
  }
  for (m=3; m-2<Mmax-1; m++)
    for (n=2; n<Nmax; n++) {
      TStbl[m-2][n] =  (TStbl[m-3][n-1] + (n-2.0*a)*TStbl[m-2][n-1]) / (n-a);
    }
}
#endif

/*
 *   stores the S table as a float
 */
#define tblSNM(N,M) S_m[M][N]
double **S_m;
double S_a;
/*
 *  dimensions
 */
int usedM, usedN;


void S_remake(double a) {
  int N, M;
  for (N=2; N<usedN; N++) {
    S_m[1][N] = log(N-a-1.0) + S_m[1][N-1];
    for (M=2; M<=usedM && M<N; M++) {
      tblSNM(N,M) = logadd(tblSNM(N-1,M-1),log(N-(M*a)-1.0)+tblSNM(N-1,M));
    }
  }    
}
  
double S_safe(int N, int T) {
  if ( N==T )
    return 0;
  return tblSNM(N,T);
}
/*
 *  fills S table
 */
void S_make(int maxN, int maxM, double a) {
  int N, M;

  usedN = maxN;
  usedM = maxM;
  
  S_m = malloc(sizeof(S_m[0])*(maxM+1));
  S_m[0] = NULL;
  for (M=1; M<=maxM; M++) {
    S_m[M] = malloc(sizeof(S_m[0][0])*maxN);
  }

  /*
   *  all values outside bounds to log(0);
   *  when N=M set to log(1)
   */
  for (N=1; N<maxN; N++) {
    if ( N<=maxM )
        tblSNM(N,N) = 0;
  }

  S_m[1][1] = 0;
  S_remake(a);
}

/*
 *   stores the S table as a float
 */
#define ftblSNM(N,M) fS_m[M][N]
float **fS_m;
float fS_a;


void fS_remake(double a) {
  int N, M;
  for (N=2; N<usedN; N++) {
    fS_m[1][N] = log(N-a-1.0) + fS_m[1][N-1];
    for (M=2; M<=usedM && M<N; M++) {
      ftblSNM(N,M) = flogadd(ftblSNM(N-1,M-1),log(N-(M*a)-1.0)+ftblSNM(N-1,M));
    }
  }    
}
  
float fS_safe(int N, int T) {
  if ( N==T )
    return 0;
  return ftblSNM(N,T);
}
/*
 *  fills S table
 */
void fS_make(float a) {
  int N, M;
  int maxN=usedN;
  int maxM=usedM;

  fS_m = malloc(sizeof(fS_m[0])*(maxM+1));
  fS_m[0] = NULL;
  for (M=1; M<=maxM; M++) {
    fS_m[M] = malloc(sizeof(fS_m[0][0])*maxN);
  }

  /*
   *  all values outside bounds to log(0);
   *  when N=M set to log(1)
   */
  for (N=1; N<maxN; N++) {
    if ( N<=maxM )
        ftblSNM(N,N) = 0;
  }

  fS_m[1][1] = 0;
  fS_remake(a);
}

double T_a = 0;
int T_Nmax=0, T_Mmax=0;
/*
 *   stores  V^{n-1}_{m-2,a} 
 */
double **Ttbl=NULL;

void T_remake(double a) {
  int n, m;
  T_a = a;
  /*
   *  m=2 case for n>=m
   */
  Ttbl[0][1] = 1.0/(1.0-a);
  for (n=2; n<T_Nmax; n++)
    Ttbl[0][n] = (1.0+(n-2*a)*Ttbl[0][n-1])/(n-a);
  /*
   *  m>2 case
   */
  Ttbl[1][2] = 1/3.0/(1.0-a);
  for (n=3; n<T_Nmax; n++) {
    for (m=3; m<=n && m<=T_Mmax; m++)
      Ttbl[m-2][n] = (1.0+(n-m*a)*Ttbl[m-2][n-1])
	/ (1.0/Ttbl[m-3][n-1]+(n-(m-1)*a));
    if ( T_Mmax>=n+1 ) 
      Ttbl[n-1][n] =  1.0 / (1.0/Ttbl[n-2][n-1]+n*(1.0-a));
  }
}

void T_make(int N, int M, double a) {
  int n;
  T_Nmax = N;
  T_Mmax = M;
  
  Ttbl = malloc((M-1)*sizeof(*Ttbl));
  Ttbl[0] = malloc(N*(M-1)*sizeof(**Ttbl));
  for (n=1; n<M-1; n++) 
    Ttbl[n] = Ttbl[n-1] + N;

  T_remake(a);
}

/*
 *   return T^n_{m,a} = S^{n}_{m,a}/S^{n}_{m-1,a} 
 *   from tables for m>1
 */
double T_V(int n, int m) {
  assert(m>1);
  assert(n>=m);
  assert(n<=T_Nmax);
  assert(m<=T_Mmax);
  return Ttbl[m-2][n-1];
}

float **fTtbl=NULL;

void fT_remake(double a) {
  int n, m;
  /*
   *  m=2 case for n>=m
   */
  fTtbl[0][1] = 1.0/(1.0-a);
  for (n=2; n<T_Nmax; n++)
    fTtbl[0][n] = (1.0+(n-2*a)*fTtbl[0][n-1])/(n-a);
  /*
   *  m>2 case
   */
  fTtbl[1][2] = 1/3.0/(1.0-a);
  for (n=3; n<T_Nmax; n++) {
    for (m=3; m<=n && m<=T_Mmax; m++)
      fTtbl[m-2][n] = (1.0+(n-m*a)*fTtbl[m-2][n-1])
	/ (1.0/fTtbl[m-3][n-1]+(n-(m-1)*a));
    if ( T_Mmax>=n+1 ) 
      fTtbl[n-1][n] =  1.0 / (1.0/fTtbl[n-2][n-1]+n*(1.0-a));
  }
}

void fT_make(int N, int M, float a) {
  int n;
  
  fTtbl = malloc((M-1)*sizeof(*fTtbl));
  fTtbl[0] = malloc(N*(M-1)*sizeof(**fTtbl));
  for (n=1; n<M-1; n++) 
    fTtbl[n] = fTtbl[n-1] + N;

  fT_remake(a);
}

/*
 *   return T^n_{m,a} = S^{n}_{m,a}/S^{n}_{m-1,a} 
 *   from tables for m>1
 */
float fT_V(int n, int m) {
  assert(m>1);
  assert(n>=m);
  assert(n<=T_Nmax);
  assert(m<=T_Mmax);
  return fTtbl[m-2][n-1];
}


int main() {
        struct timespec tstart, tend;
	int m;
	int DIM=10000+1;
	int MM=4000;
#ifdef try_TS
	if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
           exit(1);
	TS_make(DIM,MM,0.5);
	if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
           exit(1);
	printf("TS_make(%d,%d) took %lf ms\n", DIM, MM, timespecDiff(&tend, &tstart)/1.0e6);
#endif
	if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
           exit(1);
	S_make(DIM,MM,0.5);
	if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
           exit(1);
	printf("S_make(%d,%d) took %lf ms\n", DIM, MM, timespecDiff(&tend, &tstart)/1.0e6);

	fS_make(0.5);
	if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
           exit(1);
	T_make(DIM,MM,0.5);
	if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
           exit(1);
	printf("T_make(%d,%d) took %lf ms\n", DIM, MM, timespecDiff(&tend, &tstart)/1.0e6);
	fT_make(DIM,MM,0.5);
	
	for (m=10; m<MM; m+=10) 
		printf("(%d,%d):  S=%10lg  fS=%10g   T=%10lg   fT=%10g\n",
			DIM-1, m, 
			exp(S_safe(DIM-1,m)-S_safe(DIM-1,m-1)),
			exp(fS_safe(DIM-1,m)-fS_safe(DIM-1,m-1)),
			T_V(DIM-1,m),
			fT_V(DIM-1,m) );
#ifdef try_TS
	for (m=10; m<MM; m+=10) 
	  printf("(%d,%d):  S=%10lg   TS=%10lg  tbl=%10lg\n",
		 DIM-1, m, 
		 S_safe(DIM-1,m),
		 TS_safe(DIM-1,m),
		 TStbl[m-2][DIM-2] );
#endif
	return 1;
}
