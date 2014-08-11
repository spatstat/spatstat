/*
   linalg.c

   Home made linear algebra

   Yes, really

   $Revision: 1.3 $ $Date: 2011/11/20 03:50:17 $ 

*/

#include <R.h>
#include <R_ext/Utils.h>

/*
    computes the weighted sum of outer products 
    y = sum (w[j] * x[,j] %o% x[,j])
*/

void wsumouter(x, n, p, w, y) 
  double *x;    /* p by n matrix */
  int *n, *p;
  double *w;    /* weight vector, length n */
  double *y;    /* output matrix p by p, initialised to zero */
{
  int N, P;
  register int i, j, k;
  double wj, xij, wjxij, xkj;
  double *xcolj;
  N = *n; 
  P = *p;
  for(j = 0; j < N; j++) {
    R_CheckUserInterrupt();
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
  register int i, j, k;
  double xij, wjxij, xkj, vik, yj;
  double *xcolj;
  N = *n; 
  P = *p;
  for(j = 0; j < N; j++) {
    R_CheckUserInterrupt();
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
