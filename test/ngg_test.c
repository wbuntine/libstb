/*
 *   This is a simple test.
 *   (C) Wray Buntine 2011
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 *   Compile with:
 *        cc -o ngg ngg_test.c -lm -lrt -lgsl -lgslcblas 
 *   then run
 *        ./ngg
 *
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>


/* Continued fraction which occurs in evaluation
 * of Q(a,x) or Gamma(a,x).
 *
 * 1 (1-a)/x 1/x (2-a)/x 2/x (3-a)/x
 * F(a,x) = ---- ------- ----- -------- ----- -------- ...
 * 1 + 1 + 1 + 1 + 1 + 1 +
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
 *
 */
static double
gamma_inc_F_CF(const double a, const double x)
{
 const int nmax = 5000;
 const double small = gsl_pow_3 (GSL_DBL_EPSILON);

 double hn = 1.0; /* convergent */
 double Cn = 1.0 / small;
 double Dn = 1.0;
 int n;

 /* n == 1 has a_1, b_1, b_0 independent of a,x,
    so that has been done by hand */
 for ( n = 2 ; n < nmax ; n++ )
   {
     double an;
     double delta;
     
     if(GSL_IS_ODD(n))
       an = 0.5*(n-1)/x;
     else
       an = (0.5*n-a)/x;
     
     Dn = 1.0 + an * Dn;
     if ( fabs(Dn) < small )
       Dn = small;
     Cn = 1.0 + an/Cn;
     if ( fabs(Cn) < small )
       Cn = small;
     Dn = 1.0 / Dn;
     delta = Cn * Dn;
     hn *= delta;
     if(fabs(delta-1.0) < GSL_DBL_EPSILON) break;
   }
 
 if(n == nmax)
   GSL_ERROR ("error in CF for F(a,x)", GSL_EMAXITER);
 return hn;
}

static double sf_gamma_inc(double a, double x, double *expfact) {
 double cf_res = gamma_inc_F_CF(a, x);
 *expfact = (a-1.0)*log(x) - x;
 return cf_res;
}

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

long double logdiff(long double V, long double lp) {
  if ( lp>V ) {
    fprintf(stderr,"logdiff(%lf,%lf) illegal\n", 
	    (double)V, (double)lp);
    exit(1);
  }
  return lp + log(exp(V-lp)-1.0);
}

/*
 *   stores the table of values:
 *
 *       S_mtx[K][N] = \frac{\Gamma(N)}{a^{N-1}} T^{N,K}_{a,M}
 */
long double **S_mtx;
double S_a = 0.5;
double S_M = 1;

long double *G_1mnM;
/*
 *  dimensions
 */
int usedK, usedN;

/*
 *   build the vector S_mtx[K][.] for fixed K
 *   using the standard summation of Eqn (26)
 */
void S_vector(double a, double M, int K) {
  int n, i, k;
  double f2;
  long double f1;
  
  /*
   *    fill up
   *       G_1mnM[n] = Gamma(K-n/a,M) * pow(M,n/a) * exp(M)
   */
  // printf("log G(%d-n/a,M):", K);
  for (n=0; n<=usedN; n++) {
    f1 = sf_gamma_inc(K-n/a,M,&f2);
    // printf(" %lg", (double)(log(f1)+f2+log(M)*n/a));
    /*
     *  the approximate formula is
     *    (K-n/a)*log(M) - M - log(n/a-K+M+1)
     */
    G_1mnM[n] = f1*exp(f2+log(M)*n/a+M);
  }
  // printf("\n");
  S_mtx[K][1] = G_1mnM[0];

  /*
   *   now do the summation,
   *   need to rewrite to stop instability
   */
  for (n=2; n<=usedN; n++) {
    long double fact;
    /*
     *    start in the middle and work out
     */
    // #define CLEVERSUM
#define CSP 10000
#ifdef CLEVERSUM
    if ( (n % 2)==0 ) {
      long double term ;
      int n1=0, nL = (n-1)/2;
      fact = 1;
      if ( n==CSP ) printf("S_mtx[%d][%d]: ", K, n);
      S_mtx[K][n] = fact * (G_1mnM[n1]-G_1mnM[n-1-n1]);
      n1++;
      fact *= ((long double)n-n1)/n1;
      for (; n1<nL;  ) {
	term = -fact * (G_1mnM[n1]-G_1mnM[n-1-n1]);
	n1++;
	fact *= ((long double)n-n1)/n1;
	term += fact * (G_1mnM[n1]-G_1mnM[n-1-n1]);
	n1++;
	fact *= ((long double)n-n1)/n1;
	S_mtx[K][n] += term;
	if ( n==CSP ) printf("+%lg(%lg,%lg) ", (double)S_mtx[K][n], 
			    (double)fact,
			    (double)(G_1mnM[n1]-G_1mnM[n-1-n1])
			    );
      }
      if ( n1==nL ) {
	S_mtx[K][n] -= fact * (G_1mnM[n1]-G_1mnM[n-1-n1]);
	if ( n==CSP ) printf("+%lg(%lg,%lg) ", (double)S_mtx[K][n], 
			    (double)fact,
			    (double)(G_1mnM[n1]-G_1mnM[n-1-n1])
			    );
      }
      if ( n==CSP ) printf("\n");
     } else {
      long double term ;
      int n1=0, nL = (n-1)/2-1;
      fact = 1;
      S_mtx[K][n] = fact * (G_1mnM[n1]+G_1mnM[n-1-n1]-2*G_1mnM[0]);
      n1++;
      fact *= ((long double)n-n1)/n1;
      for (; n1<nL; ) {
	term = -fact * (G_1mnM[n1]+G_1mnM[n-1-n1]-2*G_1mnM[0]);
	n1++;
	fact *= ((long double)n-n1)/n1;
	term += fact * (G_1mnM[n1]+G_1mnM[n-1-n1]-2*G_1mnM[0]);
	n1++;
	fact *= ((long double)n-n1)/n1;
	S_mtx[K][n] += term;
	if ( n==CSP ) printf("+%lg(%lg) ", (double)S_mtx[K][n], (double)fact);
      }
      if ( n1==nL ) {
	S_mtx[K][n] -= fact * (G_1mnM[n1]+G_1mnM[n-1-n1]-2*G_1mnM[0]);
	if ( n==CSP ) printf("+%lg(%lg) ", (double)S_mtx[K][n], (double)fact);
	n1++;
	fact *= -((long double)n-n1)/n1;
      }
      S_mtx[K][n] -= fact * (G_1mnM[n1]-G_1mnM[0]);
      if ( n==CSP ) printf("+%lg(%lg) ", (double)S_mtx[K][n], (double)fact);
    }
    // #define TESTCLEVER
#ifdef TESTCLEVER
    {
      long double TRY = 0;
      fact = 1;
      for (i=1; i<=n; i++) {
	fact *= -1.0 * (n-i) / i;
	TRY += fact*(G_1mnM[i]-G_1mnM[0]);
      }
      printf("T^{%d,%d}*e(-M) Old = %10lg, New = %10lg\n",
	     n, K, (double)TRY, (double)S_mtx[K][n]);
    }
#endif
#else
    fact = 1;
    S_mtx[K][n] = G_1mnM[0];
    for (i=1; i<=n; i++) {
      fact *= -1.0 * (n-i) / i;
      S_mtx[K][n] += fact*G_1mnM[i];
    }
#endif
  }
  for (n=1; n<=usedN; n++) {
    /*
     *    put back "-M" since added "-M" originally
     */
    // printf("T^{%d,%d}*e(-M): %lg \n", n, K, (double)S_mtx[K][n]);
    S_mtx[K][n] = log(S_mtx[K][n]) - M;
    // printf("T^{%d,%d}: %lg \n", n, K, (double)S_mtx[K][n]);
  }
}

void S_remake(double a, double M) {
  int n, i, k;
  
  /*
   *    build base case S_mtx[1][..]
   */
  S_vector(a,M,1);
  // #define TEST1 3
#ifdef TEST1
  {
    int jj = 2;
    for ( jj=2; jj<=TEST1 ; jj++)
      S_vector(a,M,jj);
    for (n=1; n<=usedN; n++) {
      printf("T^{%d,...}: ", n);
      for (k=1; k<=TEST1 && k<=n; k++) {
	printf(" %lg", (double)S_mtx[k][n]);
      }
      printf("\n");
    }
  }
#endif

  for (n=2; n<=usedN; n++) {
    for (k=2; k<=usedK && k<=n; k++) {
#define REPORT1
#ifdef REPORT1
      printf("T^(%d,%d): T^{%d,%d}=%lg T^{%d,%d}=%lg %lg", n, k,
	     n, k-1, (double)S_mtx[k-1][n], 
	     n-1, k-1, (double)S_mtx[k-1][n-1], 
	     (double)(log((n-1)/a-(k-1.0))+S_mtx[k-1][n]) );
#endif
      /*
       *  looks different to standard recursion because
       *  of the missing a^(N-1)/\Gamma(N) term
       */
      S_mtx[k][n] = logdiff(log((n-1.0)/a) + S_mtx[k-1][n-1],
			    log((n-1.0)/a-(k-1.0))+S_mtx[k-1][n]);
#ifdef REPORT1
      printf(" -> %lg\n", (double)S_mtx[k][n]);
#endif
    }
  }    
}

#define ALTERNATE_INTEGRAL
#ifdef ALTERNATE_INTEGRAL
double smax(double a, double M, int N, int K) {
  int i;
  double s=0.5;
  double sigma;
  double logarea;
#if 0
  int STEPS=100;
  printf("dS:  ");
  for (i=0; i<STEPS; i++) {
    s = (i+0.5)/(STEPS+1);
    printf(" %lf", 
	   (N-1)/s + (1+K*a)/(1-s) - a*M/pow(1-s,1+a)
	   );
  }
  printf("\n");
  printf("S:  ");
  for (i=0; i<STEPS; i++) {
    s = (i+0.5)/(STEPS+1);
    printf(" %lf", 
	   (N-1)*log(s) - (1+K*a)*log(1-s) + M*(1-pow(1-s,-a))
	   );
  }
  printf("\n");
#endif
  /*
   *    finds the location of maximum for 3rd integral rep.
   *    i.e.,  s = u/(1+u)  for s \in [0,1)
   */
  s = 0.5;
#if 0
  printf("N-R S:  ");
#endif
  for (i=0; i<20; i++) {
    double delta = ((N-1)/s + (1+K*a)/(1-s) - a*M/pow(1-s,1+a));
    double dd = (-(N-1)/s/s + (1+K*a)/(1-s)/(1-s) - a*(1+a)*M/pow(1-s,2+a));
    if ( dd >=0 ) {
      printf("smax:  2nd deriv non neg");
      exit(1);
    }
    if ( delta<0 ) {
      delta /= -dd;
      if ( delta < -0.1 )
	delta = -0.1;
      if ( delta*1.1+s<=0 )
	s /=2;
      else
	s += delta;
    } else {
      delta /= -dd;
      if ( delta > 0.1 )
	delta = 0.1;
      if ( delta*1.1+s>=1 )
	s += (1-s)/2;
      else
	s += delta;
    }
#if 0
    printf(" %lf(%lg)", s, 1.0/sqrt(-dd));
#endif
    sigma = -dd;
  }
  /*
   *  Gaussian approx:
   *     area = prob(\hat{s}) * sqrt(2 \pi \hat\sigma^2)
   *  where
   *     \hat\sigma^2 = - \left( d^2 \log prob(s)/ ds^2 \right)^{-1}
   */
  logarea = 0.5 * log(2*M_PI/sigma) + 
    (N-1)*log(s) - (1+K*a)*log(1-s) + M*(1-pow(1-s,-a));
#if 0
  printf("Log area =~ %lg\n", logarea);
#endif
  return logarea;
}

double umax(double a, double M, int N, int K) {
  int i;
  double eu, u=0.1;
  double sigma;
  double logarea;
  double delta, dd;
  /*
   *    finds the location of maximum for 2nd integral rep.
   */
#if 0
  printf("N-R U(%d,%d):  ", N, K);
#endif
  for (i=0; i<20; i++) {
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
#if 0
    printf(" %lg(%lg,%lg)", u, delta,dd);
#endif
    sigma = dd;
  }
  /*
   *  Gaussian approx:
   *     area = prob(\hat{u}) * sqrt(2 \pi \hat\sigma^2)
   *  where
   *     \hat\sigma^2 = - \left( d^2 \log prob(u)/ du^2 \right)^{-1}
   */
  eu = exp(u);
  logarea = 0.5 * log(2*M_PI/sigma) -
    ((N-K*a)*log(1+eu) - N*u  - M * (1-pow(1+eu,a)));
#if 0
  printf("Log area =~ %lg\n", logarea);
#endif
  return logarea;
}
#endif

double S_safe(int N, int K) {
  return S_mtx[K][N];
}
/*
 *  fills S table
 */
void S_make(int maxN, int maxK, double a, double M) {
  int N, K;
  
  usedN = maxN;
  usedK = maxK;
  
  S_mtx = malloc(sizeof(S_mtx[0])*(maxK+1));
  G_1mnM = malloc(sizeof(G_1mnM[0])*(maxN+1));
  S_mtx[0] = NULL;
  for (K=1; K<=maxK; K++) {
    S_mtx[K] = malloc(sizeof(S_mtx[0][0])*(maxN+1));
  }
  S_a = a;
  S_M = M;
  S_remake(a,M);
}
 
double T_a = 0.5;
double T_M = 1;
int T_Nmax=0, T_Kmax=0;
/*
 *   stores  V^{n-1}_{k-2,a} 
 */
double **Ttbl=NULL;

void T_remake(double a, double M) {
  int n, k;
  /*
   *  base case, Ttbl[1][n]
   */
  for (n=1; n<=T_Nmax; n++)
    Ttbl[1][n] = exp(S_mtx[1][n+1]-S_mtx[1][n]);
  /*
   *  k>=2 case
   */
  for (k=1; k<T_Kmax; k++) {
    for (n=k+1; n<=T_Nmax; n++) {
      Ttbl[k+1][n] = (k + (1.0-Ttbl[k][n])*n/a)
	/ (k + (1.0/Ttbl[k][n-1]-1.0)*(n-1)/a);
    }
  }
}

void T_make(int N, int K, double a, double M) {
  int n;
  T_Nmax = N;
  T_Kmax = K;
  
  Ttbl = malloc((K+1)*sizeof(*Ttbl));
  Ttbl[0] = malloc((N+1)*(K+1)*sizeof(**Ttbl));
  for (n=1; n<=K; n++) 
    Ttbl[n] = Ttbl[n-1] + (N+1);
  
  T_a = a;
  T_M = M;
  T_remake(a, M);
}

/*
 *   return \Lambda_1 T^{n,k}_{k,a}
 */
double T_L1(int n, int k) {
  assert(k>=1);
  assert(n>=k);
  assert(n<=T_Nmax);
  assert(k<=T_Kmax);
  return Ttbl[k][n];
}
double T_L2(int n, int k) {
  assert(k>=1);
  assert(n>=2);
  assert(n<=T_Nmax);
  assert(k<=T_Kmax);
  return 1.0/Ttbl[k][n-1] + (k-n/T_a);
}

int main() {
        struct timespec tstart, tend;
	double p_a = 0.1, p_M = 10;
	int n, k;
	int MAXN=15+1;    //   maximum for N
	int MAXK=10;       //   maximum for K

	if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
           exit(1);
	S_make(MAXN,MAXK, p_a, p_M);
	if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
           exit(1);
	printf("S_make(%d,%d) took %lf ms\n", MAXN, MAXK, timespecDiff(&tend, &tstart)/1.0e6);

	if ( clock_gettime(CLOCK_MONOTONIC, &tstart)!=0 )
           exit(1);

	T_make(MAXN,MAXK, p_a, p_M);
	if ( clock_gettime(CLOCK_MONOTONIC, &tend)!=0 )
           exit(1);
	printf("T_make(%d,%d) took %lf ms\n", MAXN, MAXK, timespecDiff(&tend, &tstart)/1.0e6);
	
	for (n=1; n<MAXN; n+=1) 
	  for (k=1; k<MAXK; k+=1) 
	    if ( k<=n ) {
	      printf("(%d,%d):  T=%10lg  G=%10lg L1(T)=%10lg   L1=%10lg\n",
		     n, k, 
		     (double)S_mtx[k][n],
		     umax(p_a,p_M,n,k)+k*log(p_M) - p_M, 
		     exp((double)(S_mtx[k][n+1]-S_mtx[k][n])),
		     T_L1(n,k) );
	    }
	return 0;
}
