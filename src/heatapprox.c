/* 
   Approximation to heat kernel 

   Copyright (c) Greg McSwiggan and Adrian Baddeley 2017-2020

   $Revision: 1.1 $ $Date: 2020/04/05 03:38:35 $

*/

#include <R.h>
#include <Rmath.h>

void heatApprox(n, a, x, y, s, degl, degr, m, z)
     int    *n; /* number of calculations (length of each vector) */
     double *a; /* rod length */
     double *x; /* source position */
     double *y; /* query position */
     double *s; /* bandwidth */
     int    *degl; /* vertex degree of left endpoint */
     int    *degr; /* vertex degree of right endpoint */
     int    *m; /* number of terms in sum */
     double *z; /* result */
{
  register int i, k, N, M, dL, dR;
  register double Z, A, twoA, Bk, X, Y, sigma, fL, fR, fLfR, cc;

  N = *n;
  M = *m;

  for(i = 0; i < N; i++) {
    sigma = s[i];
    A = a[i];
    if(A <= 0.0 || sigma <= 0.0) {
      /* trap bad data */
      z[i] = 0.0;
    } else {
      /* do calculation */
      X = x[i];
      Y = y[i];
      dL = degl[i];
      dR = degr[i];
      fL = (2.0/dL - 1.0);
      fR = (2.0/dR - 1.0);
      fLfR = fL * fR;
      twoA = 2.0 * A;
      Z = dnorm(Y, X, sigma, (int) 0);
      cc = 1.0;
      for(k = 1; k <= M; k++) {
        Bk = k * twoA;
        Z += cc * (fR   * dnorm( Bk - Y, X, sigma, (int) 0)  + 
		   fLfR * dnorm( Bk + Y, X, sigma, (int) 0) + 
		   fL   * dnorm(-Bk + Y, X, sigma, (int) 0) + 
		   fLfR * dnorm(-Bk - Y, X, sigma, (int) 0));
	cc *= fLfR;
      }
      z[i] = Z;
    }
  }
}

  

  
  
  
