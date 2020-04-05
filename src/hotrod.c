/* 
   hotrod.c

   Heat kernel on a one-dimensional rod
   with either insulated ends or absorbing ends

   Copyright (c) Greg McSwiggan and Adrian Baddeley 2017-2020

   $Revision: 1.1 $ $Date: 2020/04/05 03:39:01 $

*/

#include <R.h>
#include <Rmath.h>

void hotrodInsul(n, a, x, y, s, m, z)
     int    *n; /* number of calculations (length of each vector) */
     double *a; /* rod length */
     double *x; /* source position */
     double *y; /* query position */
     double *s; /* bandwidth */
     int    *m; /* number of terms in sum */
     double *z; /* result */
{
  register int i, k, N, M;
  register double Z, A, twoA, Bk, X, Y, sigma;

  N = *n;
  M = *m;

  for(i = 0; i < N; i++) {
    sigma = s[i];
    A = a[i];
    if(A <= 0.0 || sigma <= 0.0) {
      /* trap bad data */
      z[i] = 0.0;
    } else if(sigma > 20.0 * A) {
      /* uniform density */
      z[i] = 1/A;
    } else {
      /* do calculation */
      X = x[i];
      Y = y[i];
      Z = 0.0;
      twoA = 2.0 * A;
      for(k = -M; k <= M; k++) {
        Bk = k * twoA;
        Z += dnorm( Bk + Y, X, sigma, (int) 0);
        Z += dnorm( Bk - Y, X, sigma, (int) 0);
      }
      z[i] = Z;
    }
  }
}

void hotrodAbsorb(n, a, x, y, s, m, z)
     int    *n; /* number of calculations (length of each vector) */
     double *a; /* rod length */
     double *x; /* source position */
     double *y; /* query position */
     double *s; /* bandwidth */
     int    *m; /* number of terms in sum */
     double *z; /* result */
{
  register int i, k, N, M;
  register double Z, A, X, Y, sigma, pionL, piXonL, piYonL, fac;

  N = *n;
  M = *m;

  for(i = 0; i < N; i++) {
    sigma = s[i];
    A = a[i];
    if(A <= 0.0 || sigma <= 0.0 || sigma > 20.0 * A) {
      z[i] = 0.0;
    } else {
      /* do calculation */
      pionL = M_PI/A;
      fac = pionL * pionL * sigma * sigma/2.0;
      X = x[i];
      Y = y[i];
      piXonL = pionL * X;
      piYonL = pionL * Y;
      Z = 0.0;
      for(k = 1; k <= M; k++) {
        Z += exp(-fac * k * k) * sin(k * piXonL) * sin(k * piYonL);
      }
      Z *= 2.0/A;
      z[i] = Z;
    }
  }
}

