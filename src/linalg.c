/*
   linalg.c

   Home made linear algebra

   Yes, really

   $Revision: 1.9 $ $Date: 2013/09/25 06:07:24 $ 

   Csumouter
   Cwsumouter
   Cquadform
   Csumsymouter
   Cwsumsymouter
*/

#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/* ............... matrices ..............................*/

/*
    Csumouter
    computes the sum of outer products of columns of x
    y = sum[j] (x[,j] %o% x[,j])
*/

void Csumouter(x, n, p, y) 
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
    Cwsumouter
    computes the weighted sum of outer products of columns of x
    y = sum[j] (w[j] * x[,j] %o% x[,j])
*/

void Cwsumouter(x, n, p, w, y) 
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

void Cquadform(x, n, p, v, y) 
  double *x;    /* p by n matrix */
  int *n, *p;
  double *v;    /* p by p matrix */
  double *y;    /* output vector, length n */
{
  int N, P;
  register int i, j, k, maxchunk;
  register double xij, xkj, vik, yj;
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

/*
    computes the bilinear form values
    z[j] = x[,j] %*% v %*% t(y[,j])
*/

void Cbiform(x, y, n, p, v, z) 
     double *x, *y;    /* p by n matrices */
     int *n, *p;
     double *v;    /* p by p matrix */
     double *z;    /* output vector, length n */
{
  int N, P;
  register int i, j, k, maxchunk;
  register double xij, vik, ykj, zj;
  register double *xcolj, *ycolj;
  N = *n; 
  P = *p;
  OUTERCHUNKLOOP(j, N, maxchunk, 2048) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, N, maxchunk, 2048) {
      xcolj = x + j * P;
      ycolj = y + j * P;
      zj = 0;
      for(i = 0; i < P; i++) {
	xij = xcolj[i];
	for(k = 0; k < P; k++) {
	  ykj = ycolj[k];
	  vik = v[k * P + i];
	  zj += xij * vik * ykj;
	}
      }
      z[j] = zj;
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

#define FNAME Csumsymouter
#include "sumsymouter.h"
#undef FNAME

/*
  wsumsymouter
  computes the weighted sum of outer products 
  w[i,j] * (x[,i,j] %o% x[,j,i])  over all pairs i, j
*/

#define FNAME Cwsumsymouter
#define WEIGHTED
#include "sumsymouter.h"
#undef FNAME
#undef WEIGHTED
