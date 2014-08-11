/*
   linalg.c

   Home made linear algebra

   Yes, really

   $Revision: 1.6 $ $Date: 2013/04/18 08:22:33 $ 

   sumouter
   wsumouter
   quadform
   sumsymouter
   wsumsymouter
*/

#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/* ............... matrices ..............................*/

/*
    sumouter
    computes the sum of outer products of columns of x
    y = sum[j] (x[,j] %o% x[,j])
*/

void sumouter(x, n, p, y) 
  double *x;    /* p by n matrix */
  int *n, *p;
  double *y;    /* output matrix p by p, initialised to zero */
{
  int N, P;
  register int i, j, k, maxchunk;
  register double xij, xkj;
  register double *xcolj;
  N = *n; 
  P = *p;
  OUTERCHUNKLOOP(j, N, maxchunk, 2048) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, N, maxchunk, 2048) {
      xcolj = x + j * P;
      for(i = 0; i < P; i++) {
	xij = xcolj[i];
	for(k = 0; k < P; k++) {
	  xkj = xcolj[k];
	  y[k * P + i] += xij * xkj;
	}
      }
    }
  }
}

/*
    wsumouter
    computes the weighted sum of outer products of columns of x
    y = sum[j] (w[j] * x[,j] %o% x[,j])
*/

void wsumouter(x, n, p, w, y) 
  double *x;    /* p by n matrix */
  int *n, *p;
  double *w;    /* weight vector, length n */
  double *y;    /* output matrix p by p, initialised to zero */
{
  int N, P;
  register int i, j, k, maxchunk;
  register double wj, xij, wjxij, xkj;
  register double *xcolj;
  N = *n; 
  P = *p;
  OUTERCHUNKLOOP(j, N, maxchunk, 2048) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, N, maxchunk, 2048) {
      wj = w[j];
      xcolj = x + j * P;
      for(i = 0; i < P; i++) {
	xij = xcolj[i];
	wjxij = wj * xij;
	for(k = 0; k < P; k++) {
	  xkj = xcolj[k];
	  y[k * P + i] += wjxij * xkj;
	}
      }
    }
  }
}

/*
    computes the quadratic form values
    y[j] = x[,j] %*% v %*% t(x[,j])
*/

void quadform(x, n, p, v, y) 
  double *x;    /* p by n matrix */
  int *n, *p;
  double *v;    /* p by p matrix */
  double *y;    /* output vector, length n */
{
  int N, P;
  register int i, j, k, maxchunk;
  register double xij, wjxij, xkj, vik, yj;
  register double *xcolj;
  N = *n; 
  P = *p;
  OUTERCHUNKLOOP(j, N, maxchunk, 2048) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, N, maxchunk, 2048) {
      xcolj = x + j * P;
      yj = 0;
      for(i = 0; i < P; i++) {
	xij = xcolj[i];
	for(k = 0; k < P; k++) {
	  xkj = xcolj[k];
	  vik = v[k * P + i];
	  yj += xij * vik * xkj;
	}
      }
      y[j] = yj;
    }
  }
}

/* ............... 3D arrays ...................... */

#undef FNAME
#undef WEIGHTED

/*
  sumsymouter
  computes the sum of outer products 
  x[,i,j] %o% x[,j,i]  over all pairs i, j
*/

#define FNAME sumsymouter
#include "sumsymouter.h"
#undef FNAME

/*
  wsumsymouter
  computes the weighted sum of outer products 
  w[i,j] * (x[,i,j] %o% x[,j,i])  over all pairs i, j
*/

#define FNAME wsumsymouter
#define WEIGHTED
#include "sumsymouter.h"
#undef FNAME
#undef WEIGHTED
