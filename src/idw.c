/*
  idw.c

  Inverse-distance weighted smoothing

  $Revision: 1.8 $ $Date: 2013/05/27 02:09:10 $

*/

#include <Rmath.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

#define MAT(X,I,J,NROW) (X)[(J) + (NROW) * (I)]

/*  inverse-distance smoothing from data points onto pixel grid */

void Cidw(x, y, v, n, xstart, xstep, nx, ystart, ystep, ny, power, num, den, rat)
     double *x, *y, *v;           /* data points and values */
     int *n;
     double *xstart, *xstep, *ystart, *ystep;   /* pixel grid */
     int *nx, *ny;
     double *power;                   /* exponent for IDW */
     double *num, *den, *rat;     /* output arrays - assumed initialised 0 */
{
  int N, i, Nx, Ny, ix, iy;
  double xg, yg, x0, dx, y0, dy, pon2, d2, w;
  
  N  = *n;
  Nx = *nx;
  Ny = *ny;
  x0 = *xstart;
  y0 = *ystart;
  dx = *xstep;
  dy = *ystep;

  pon2 = (*power)/2.0;

  if(pon2 == 1.0) {
    /* slightly faster code when power=2 */
    for(ix = 0, xg=x0; ix < Nx; ix++, xg+=dx) {
      if(ix % 256 == 0) R_CheckUserInterrupt();
      for(iy = 0, yg=y0; iy < Ny; iy++, yg+=dy) {
	/* loop over data points, accumulating numerator and denominator */
	for(i = 0; i < N; i++) {
	  d2 = (xg - x[i]) * (xg - x[i]) + (yg - y[i]) * (yg - y[i]);
	  w = 1.0/d2;
	  MAT(num, ix, iy, Ny) += w * v[i];
	  MAT(den, ix, iy, Ny) += w;
	}
	/* compute ratio */
	MAT(rat, ix, iy, Ny) = MAT(num, ix, iy, Ny)/MAT(den, ix, iy, Ny);
      }
    }
  } else {
    /* general case */
    for(ix = 0, xg=x0; ix < Nx; ix++, xg+=dx) {
      if(ix % 256 == 0) R_CheckUserInterrupt();
      for(iy = 0, yg=y0; iy < Ny; iy++, yg+=dy) {
	/* loop over data points, accumulating numerator and denominator */
	for(i = 0; i < N; i++) {
	  d2 = (xg - x[i]) * (xg - x[i]) + (yg - y[i]) * (yg - y[i]);
	  w = 1.0/pow(d2, pon2);
	  MAT(num, ix, iy, Ny) += w * v[i];
	  MAT(den, ix, iy, Ny) += w;
	}
	/* compute ratio */
	MAT(rat, ix, iy, Ny) = MAT(num, ix, iy, Ny)/MAT(den, ix, iy, Ny);
      }
    }
  }
}

/* Leave-one-out IDW at data points only */

void idwloo(x, y, v, n, power, num, den, rat)
     double *x, *y, *v;           /* data points and values */
     int *n;
     double *power;                   /* exponent for IDW */
     double *num, *den, *rat;     /* output vectors - assumed initialised 0 */
{
  int N, i, j, maxchunk;
  double xi, yi, d2, w, pon2;
  
  N  = *n;
  pon2 = (*power)/2.0;

  if(pon2 == 1.0) {
    /* slightly faster code when power=2 */
    OUTERCHUNKLOOP(i, N, maxchunk, 16384) {
      R_CheckUserInterrupt();
      INNERCHUNKLOOP(i, N, maxchunk, 16384) {
	xi = x[i];
	yi = y[i];
	if(i > 0) {
	  for(j = 0; j < i; j++) {
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/d2;
	    num[i] += w * v[j];
	    den[i] += w;
	  }
	}
	if(i < N-1) {
	  for(j = i+1; j < N; j++) {
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/d2;
	    num[i] += w * v[j];
	    den[i] += w;
	  }
	}
	/* compute ratio */
	rat[i] = num[i]/den[i];
      }
    }
  } else {
    /* general case */
    OUTERCHUNKLOOP(i, N, maxchunk, 16384) {
      R_CheckUserInterrupt();
      INNERCHUNKLOOP(i, N, maxchunk, 16384) {
	xi = x[i];
	yi = y[i];
	if(i > 0) {
	  for(j = 0; j < i; j++) {
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/pow(d2, pon2);
	    num[i] += w * v[j];
	    den[i] += w;
	  }
	}
	if(i < N-1) {
	  for(j = i+1; j < N; j++) {
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/pow(d2, pon2);
	    num[i] += w * v[j];
	    den[i] += w;
	  }
	}
	/* compute ratio */
	rat[i] = num[i]/den[i];
      }
    }
  }
}



