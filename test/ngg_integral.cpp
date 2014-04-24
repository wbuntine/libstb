/*
 *  Tester for playing with normalising NGGs
 *  Author Changyou Chen (cchangyou@gmail.com)
 *  some mods by Wray Buntine
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "DEIntegrator.h"

/*
 *   key parameters
 */
int RUNS, STEP;   
int K;
double N, a, M;
double relerror;
double interror;

/***********************************************************
 *    service routines for hmax()
 * 
 *  take integral in $u$ after (23) in report
 *  use  x = u/(1+u)
 *  Then compute    T^{N,K}_{a,M} \Gamma(N)/a^{N-1}
 ***********************************************************/
double hval_min, hval_max;
static double h(double x) {
  return - (N-1)*log(x) + (K*a+1)*log(1-x) + M/pow(1-x,a);
}
static double dh(double x) {
  return - (N-1)/x - (K*a+1)/(1-x) + M*a/pow(1-x,1+a);
}
static double ddh(double x) {
  return - (K*a+1)/(1-x)/(1-x) + (N-1)/x/x + M*a*(1+a)/pow(1-x,2+a);
}

static double hmax(double *initu) {
  int i;
  double eu, olddelta;
  double sigma;
  double logarea;
  double delta, dd;
  /*
   *    finds the location of maximum for [0,1] integral rep.
   */
  eu = *initu;
#if 0
  printf("N-R U(%lf,%d) for h(%lg) -> %lg\n", N, K, eu, h(eu));
#endif
  olddelta = 0.1;
  for (i=0; i<100; i++) {
    delta = dh(eu);
    dd = ddh(eu);
    if ( dd <=0 ) {
      int j;
#if 0
      for (eu=0.1; eu<1; eu+=0.1)
	printf("  hmax:  2nd deriv = %lg, 1st deriv = %lg\n", ddh(eu), dh(eu));
      exit(1);
#endif
      //  direction
      olddelta = fabs(olddelta);
      if ( delta<0 ) 
	olddelta *= -1;
      delta = olddelta;
      //   bounds
      if ( eu-delta<0 )
	delta = eu*0.7;
      if ( eu-delta>1 )
	delta = -(1-eu)*0.7;
      j = 0;
#if 0
      printf("  loop 0:  %lg - %lg -> %lg\n", eu, delta, h(eu-delta));
#endif
      while ( h(eu-delta)>= h(eu) ) {
	delta /= 2;
#if 0
	printf("  loop %d:  %lg - %lg -> %lg\n", j+1, eu, delta, h(eu-delta));
#endif
	if ( j++>10 ) {
	  printf("hmax looped:  2nd deriv = %lg, 1st deriv = %lg\n", ddh(eu), dh(eu));
	  exit(1);
	}
      }
    } else {
      delta /= dd;
      if ( eu-delta<0 )
	delta = eu*0.7;
      if ( eu-delta>1 )
	delta = -(1-eu)*0.7;
#if 0
      printf("  NR:  %lg - %lg -> %lg\n", eu, delta, h(eu-delta));
#endif
    } 
    eu -= delta;
    if ( fabs(delta/eu)<relerror )
      break;
#if 0
    printf(" %lg(%lg,%lg)", eu, delta,dd);
#endif
    olddelta = delta*0.9;
  }
  /*
   *  Gaussian approx:
   *     area = prob(\hat{u}) * sqrt(2 \pi \hat\sigma^2)
   *  where
   *     \hat\sigma^2 = - \left( d^2 \log prob(u)/ du^2 \right)^{-1}
   */
  *initu = eu;
  sigma = ddh(eu);
  logarea = 0.5 * log(2*M_PI/sigma) - h(eu) + K*log(M) + log(a);	
  hval_min = (N>1) ? eu*pow(interror,1.0/(N-1)) : 0;
  hval_max = 1.0 - pow(pow(1-eu,-a)-log(interror)/M,-1.0/a);
#if 0
  printf("Log area (%0lf,%lg,%lg) =~ %lg (%lg,%lg)\n", 
	 N, eu, 1/sqrt(sigma), logarea,
	 hval_min, hval_max);
#endif
  return logarea;
}

static double hmax2(double *initu) {
  int i;
  double u;
  double sigma, logarea;
  /*
   *    finds the location of maximum for [0,1] integral rep.
   */
#if 0
  printf("FP U(%lf,%d):  ", N, K);
#endif
  u = *initu;
  for (i=0; i<10; i++) {
    // double ud = 1 - pow((1+K*a)/M/a + (N-1)*(1-u)/M/a/u, -1.0/a);
    double ud = 1/(1+(M*a/pow(1-u,a)-1-K*a)/(N-1));
    if ( fabs((ud-u)/u)<relerror )
      break;
#if 1
    printf(" %lg", ud);
#endif
    u = ud;
  }
  /*
   *  Gaussian approx:
   *     area = prob(\hat{u}) * sqrt(2 \pi \hat\sigma^2)
   *  where
   *     \hat\sigma^2 = - \left( d^2 \log prob(u)/ du^2 \right)^{-1}
   */
  *initu = u;
  sigma = ddh(u);
  logarea = 0.5 * log(2*M_PI/sigma) - h(u) + K*log(M);	
  hval_min = (N>1) ? u*pow(interror,1.0/(N-1)) : 0;
  hval_max = 1.0 - pow(pow(1-u,-a)-log(interror)/M,-1.0/a);
#if 1
  printf("Log area (%0lf,%lg,%lg) =~ %lg (%lg,%lg)\n", 
	 N, u, 1/sqrt(sigma), 
	 logarea,
	 hval_min, hval_max);
#endif
  return logarea;
}

/*********************************************************
 *    service routines for gmax()
 *
 *    Take Eq (22) in report and set 
 *       x = t-M
 *    Then compute    T^{N,K}_{a,M} \Gamma(N)/a^{N-1}
 *
 ********************************************************/
//   double pxm = pow(1+x/M,1/a);
static double g(double x, double pxm) {
  return x - (K-1)*log(1+x/M) - (N-1)*log(1-1/pxm);
}
static double dg(double x, double pxm) {
  return 1 - ((K-1) + (N-1)/a/(pxm-1))/(1+x/M)/M;
}
static double ddg(double x, double pxm) {
  return ( (K-1) + (N-1)/a/(pxm-1) + (N-1)*pxm/a/a/(pxm-1)/(pxm-1) )
   /(1+x/M)/M/(1+x/M)/M;
}

static double gmax(double *initu) {
  int i;
  double pxm, eu, u;
  double sigma;
  double logarea;
  double delta, dd;
  /*
   *    finds the location of maximum for 2nd integral rep.
   */
#if 0
  printf("N-R U(%lf,%d):  ", N, K);
#endif
  u = *initu;
  for (i=0; i<100; i++) {
    eu = exp(u);
    pxm = pow(1+eu/M,1/a);
    delta = eu*dg(eu,pxm) - 1;
    dd = eu*eu*ddg(eu,pxm) + eu*dg(eu,pxm);
    if ( dd <=0 ) {
      printf("gmax:  2nd deriv non neg");
      exit(1);
    }
    delta /= dd;
    if ( delta/u < -2 )
      delta = -2*u;
    else if ( delta/u > 2 )
      delta = 2*u; 
    u -= delta;
    if ( fabs(1-exp(delta))<relerror )
      break;
#if 0
    printf(" %lg(%lg,%lg)", u, delta,dd);
#endif
  }
  /*
   *  Gaussian approx:
   *     area = prob(\hat{u}) * sqrt(2 \pi \hat\sigma^2)
   *  where
   *     \hat\sigma^2 = - \left( d^2 \log prob(u)/ du^2 \right)^{-1}
   */
  eu = exp(u);
  *initu = u;
  pxm = pow(1+eu/M,1/a);
  sigma = eu*eu*ddg(eu,pxm) + eu*dg(eu,pxm);
  logarea = 0.5 * log(2*M_PI/sigma) - g(eu,pxm) + u 
    + ((K-1)*log(M) - M);				
#if 0
  printf("Log area (%0lf,%lg,%lg) =~ %lg\n", N, u, 1/sqrt(sigma), logarea);
#endif
  return logarea;
}


static double umax(double *initu) {
  int i;
  double eu, u;
  double sigma;
  double logarea;
  double delta, dd;
  /*
   *    finds the location of maximum for 2nd integral rep.
   */
#if 0
  printf("N-R U(%d,%d):  ", N, K);
#endif
  u = *initu;
  for (i=0; i<100; i++) {
    eu = exp(u);
    delta = (N-K*a)*eu/(1+eu) - N + a*M * eu / pow(1+eu,1-a);
    dd = (N-K*a)*eu/(1+eu)/(1+eu) + a*a*M * eu / pow(1+eu,1-a)
      + a*(1-a)*M * eu / pow(1+eu,2-a);
    if ( dd <=0 ) {
      printf("umax:  2nd deriv non neg");
      exit(1);
    }
    delta /= dd;
    if ( delta/u < -2 )
      delta = -2*u;
    else if ( delta/u > 2 )
      delta = 2*u; 
    u -= delta;
    if ( fabs(1-exp(delta))<relerror )
      break;
#if 0
    printf(" %lg(%lg,%lg)", u, delta,dd);
#endif
  }
  /*
   *  Gaussian approx:
   *     area = prob(\hat{u}) * sqrt(2 \pi \hat\sigma^2)
   *  where
   *     \hat\sigma^2 = - \left( d^2 \log prob(u)/ du^2 \right)^{-1}
   */
  eu = exp(u);
  *initu = u;
  sigma = (N-K*a)*eu/(1+eu)/(1+eu) + a*a*M * eu / pow(1+eu,1-a)
    + a*(1-a)*M * eu / pow(1+eu,2-a);
  logarea = 0.5 * log(2*M_PI/sigma) -
    ((N-K*a)*log(1+eu) - N*u  - M * (1-pow(1+eu,a)))
    + K*log(M) - M + log(a);
#if 0
  printf("Log area (%0lf,%lg,%lg) =~ %lg\n", N, u, 1/sqrt(sigma), logarea);
#endif
  return logarea;
}


static int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p) {
  return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
           ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}

// integration domain transformed to [0 -- 1 / M]
class Integral
{
public:
  double operator()(double x) const
  {
    return exp((N - 1) * log(1 - pow(M * x, 1 / a)) - (K + 1) * log(x) - 1 / x);
  }
};

// integration domain [M -- \infty]
class Integral1
{
public:
  double operator()(double x) const
  {
    return exp(-x + (K - 1)*log(x) + (N - 1)*log(1 - pow(M / x, 1 / a)) );
  }
};

// integration domain [vmin,vmax]
class Integral2
{
public:
	double operator()(double x) const
	{
	  return exp(-h(x) + (K-1)*log(M));
	}
};

double *est;
double *est1;
double *approx;
double *approx2;
int  *Nval;


int main(int argc, char* argv[])
{
        struct timespec tstart, tend;
	int c;
	int i;
 	Integral f;
	Integral1 f1;
	Integral2 f2;
	int evaluations;
	double errorEstimate;
	double integral, integral1, integral2;
	int diff = 0;
	int kmax = 0;
	double saveu;

	/*
	 *  set defaults
	 */
	a = 0.1;
	K = 2;
	M = 10;
	relerror=1e-6;
	interror=1e-10;
	RUNS = 100;
	STEP = 1;
	
	while ( (c=getopt(argc, argv,"da:K:kM:N:S:r:"))>=0 ) {
	  switch ( c ) {
	  case 'd':
	    diff = 1;
	    break;
	  case 'k':
	    kmax = 1;
	    break;
	  case 'r':
	    if ( !optarg || sscanf(optarg,"%lf",&relerror)!=1 ) {
	      fprintf(stderr, "Need a valid 'r' argument\n");
	      exit(1);
	    }
	    break;
	  case 'M':
	    if ( !optarg || sscanf(optarg,"%lf",&M)!=1 ) {
	      fprintf(stderr, "Need a valid 'M' argument\n");
	      exit(1);
	    }
	    break;
	  case 'a':
	    if ( !optarg || sscanf(optarg,"%lf",&a)!=1 ) {
	      fprintf(stderr, "Need a valid 'a' argument\n");
	      exit(1);
	    }
	    break;
	  case 'K':
	    if ( !optarg || sscanf(optarg,"%d",&K)!=1 ) {
	      fprintf(stderr, "Need a valid 'K' argument\n");
	      exit(1);
	    }
	    break;
	  case 'S':
	    if ( !optarg || sscanf(optarg,"%d",&STEP)!=1 ) {
	      fprintf(stderr, "Need a valid 'S' argument\n");
	      exit(1);
	    }
	    break;
	  case 'N':
	    if ( !optarg || sscanf(optarg,"%d",&RUNS)!=1 ) {
	      fprintf(stderr, "Need a valid 'N' argument\n");
	      exit(1);
	    }
	    break;
	  default:
	    fprintf(stderr, "Bad options:\n");
	    fprintf(stderr, "  -a -M -k -K -N\n");
	    fprintf(stderr, "  -r (relative error)\n");
	    fprintf(stderr, "  -S (set size)\n");
	    exit(1);
	  }
	}

	est = (double*)malloc((RUNS/STEP+2)*sizeof(double));
	est1 = (double*)malloc((RUNS/STEP+2)*sizeof(double));
	approx = (double*)malloc((RUNS/STEP+2)*sizeof(double));
	approx2 = (double*)malloc((RUNS/STEP+2)*sizeof(double));
	Nval = (int*)malloc((RUNS/STEP+2)*sizeof(int));

	fprintf(stderr, "Dimensions: a=%lf, M=%lf, K=%d, N=%d...%d, relerror=%lg\n",
		a, M, K, K, RUNS, relerror);
        if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
	  exit(1);	
	for(i=0, N = 1; N < RUNS; N+=STEP){
	  if ( kmax ) K = N;
	  Nval[i] = N;
	  est[i] = log(DEIntegrator<Integral>::Integrate(f, 0, 1 / M, 1e-40, evaluations, errorEstimate));
	  i++;
	}
        if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
	  exit(1);
	fprintf(stderr, "Time for integrating f for N=%d, +%d, %d is %lf ms\n",
		1, STEP, RUNS, timespecDiff(&tend, &tstart)/1.0e6);

        if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
	  exit(1);	
	saveu = 0.5;
	for(i=0, N = 1; N < RUNS; N+=STEP) {
	  if ( kmax ) K = N;
	  hmax(&saveu);
	  est1[i] = 
	    // log(DEIntegrator<Integral1>::Integrate(f1, M, 1000000, 1e-40, evaluations, errorEstimate));
	    log(DEIntegrator<Integral2>::Integrate(f2, hval_min, hval_max, 1e-40, evaluations, errorEstimate));
	  i++;
	}
        if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
	  exit(1);
	fprintf(stderr, "Time for integrating f1 for N=%d, +%d, %d is %lf ms\n",
		1, STEP, RUNS, timespecDiff(&tend, &tstart)/1.0e6);

        if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
	  exit(1);	
	saveu = 0.5;
	for(i=0, N = 1; N < RUNS; N+=STEP) {
	  if ( kmax ) K = N;
	  approx2[i] = hmax(&saveu);
	  i++;
	}
        if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
	  exit(1);
	fprintf(stderr, "Time for approx hmax for N=%d, +%d, %d is %lf ms\n",
		1, STEP, RUNS, timespecDiff(&tend, &tstart)/1.0e6);
 
        if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
	  exit(1);
	saveu = log(0.5);
	for (i=0, N = 1; N < RUNS; N+=STEP) {
	  if ( kmax ) K = N;
	  approx[i] = umax(&saveu);
	  i++;
	}
        if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
	  exit(1);
	fprintf(stderr, "Time for approx umax for N=%d, +%d, %d is %lf ms\n",
		1, STEP, RUNS, timespecDiff(&tend, &tstart)/1.0e6);
 
	for(i=0, N = 1; N < RUNS; N+=STEP, i++){
	  if ( !kmax && N<K ) 
	    continue;
	  if ( diff ) {
	    if ( N==1 )
	      continue;
	    printf("%d: %lg %lg %lg %lg\n", 
		   Nval[i], (est[i]-est[i-1]), (est1[i]-est1[i-1]), 
		   (approx[i]-approx[i-1]), (approx2[i]-approx2[i-1]));
	  } else 
	    printf("%d: %lg %lg %lg %lg\n", 
		   Nval[i], est[i], est1[i], approx[i], approx2[i]);
	}
	free(est);
	free(est1);
	free(approx);
	free(approx2);
	return 0;
}
