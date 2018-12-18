/*
  idw.c

  Inverse-distance weighted smoothing

  $Revision: 1.11 $ $Date: 2018/12/18 02:43:11 $

  Cidw    inverse distance smoothing from data points onto pixel grid
  idwloo  leave-one-out estimate at data points

  Cidw2   Cidw with variance estimate
  idwloo2 idwloo with variance estimate

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

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
  double xg, yg, x0, dx, y0, dy, pon2, d2, w, sumw, sumwv;
  
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
	sumwv = sumw = 0.0;
	/* loop over data points, accumulating numerator and denominator */
	for(i = 0; i < N; i++) {
	  d2 = (xg - x[i]) * (xg - x[i]) + (yg - y[i]) * (yg - y[i]);
	  w = 1.0/d2;
	  sumwv += w * v[i];
	  sumw  += w;
	}
	/* compute ratio */
	MAT(num, ix, iy, Ny) = sumwv;
	MAT(den, ix, iy, Ny) = sumw;
	MAT(rat, ix, iy, Ny) = sumwv/sumw;
      }
    }
  } else {
    /* general case */
    for(ix = 0, xg=x0; ix < Nx; ix++, xg+=dx) {
      if(ix % 256 == 0) R_CheckUserInterrupt();
      for(iy = 0, yg=y0; iy < Ny; iy++, yg+=dy) {
	sumwv = sumw = 0.0;
	/* loop over data points, accumulating numerator and denominator */
	for(i = 0; i < N; i++) {
	  d2 = (xg - x[i]) * (xg - x[i]) + (yg - y[i]) * (yg - y[i]);
	  w = 1.0/pow(d2, pon2);
	  sumwv += w * v[i];
	  sumw  += w;
	}
	/* compute ratio */
	MAT(num, ix, iy, Ny) = sumwv;
	MAT(den, ix, iy, Ny) = sumw;
	MAT(rat, ix, iy, Ny) = sumwv/sumw;
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
  double xi, yi, d2, w, pon2, sumw, sumwv;
  
  N  = *n;
  pon2 = (*power)/2.0;

  if(pon2 == 1.0) {
    /* slightly faster code when power=2 */
    OUTERCHUNKLOOP(i, N, maxchunk, 16384) {
      R_CheckUserInterrupt();
      INNERCHUNKLOOP(i, N, maxchunk, 16384) {
	xi = x[i];
	yi = y[i];
	sumwv = sumw = 0.0;
	if(i > 0) {
	  for(j = 0; j < i; j++) {
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/d2;
	    sumwv += w * v[j];
	    sumw  += w;
	  }
	}
	if(i < N-1) {
	  for(j = i+1; j < N; j++) {
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/d2;
	    sumwv += w * v[j];
	    sumw  += w;
	  }
	}
	/* compute ratio */
	num[i] = sumwv;
	den[i] = sumw;
	rat[i] = sumwv/sumw;
      }
    }
  } else {
    /* general case */
    OUTERCHUNKLOOP(i, N, maxchunk, 16384) {
      R_CheckUserInterrupt();
      INNERCHUNKLOOP(i, N, maxchunk, 16384) {
	xi = x[i];
	yi = y[i];
	sumwv = sumw = 0.0;
	if(i > 0) {
	  for(j = 0; j < i; j++) {
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/pow(d2, pon2);
	    sumwv += w * v[j];
	    sumw  += w;
	  }
	}
	if(i < N-1) {
	  for(j = i+1; j < N; j++) {
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/pow(d2, pon2);
	    sumwv += w * v[j];
	    sumw  += w;
	  }
	}
	/* compute ratio */
	num[i] = sumwv;
	den[i] = sumw;
	rat[i] = sumwv/sumw;
      }
    }
  }
}

/* ----------------------------------------------------
   VERSIONS WITH VARIANCE CALCULATION
   --------------------------------------------------- */

/*  inverse-distance smoothing from data points onto pixel grid */

void Cidw2(x, y, v, n, xstart, xstep, nx, ystart, ystep, ny, power,
	   num, den, rat, mtwo, wtwo)
     double *x, *y, *v;           /* data points and values */
     int *n;
     double *xstart, *xstep, *ystart, *ystep;   /* pixel grid */
     int *nx, *ny;
     double *power;                   /* exponent for IDW */
     double *num, *den, *rat;    /* output arrays - assumed initialised 0 */
     double *mtwo, *wtwo;    /* output arrays - assumed initialised 0 */
{
  int N, i, Nx, Ny, ix, iy;
  double xg, yg, x0, dx, y0, dy, pon2, d2, w, vi,
    sumw, sumwv, sumw2, runmean, m2, delta, epsilon;
  
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
	sumwv = sumw = sumw2 = m2 = runmean = 0.0;
	/* loop over data points, accumulating numerator and denominator */
	for(i = 0; i < N; i++) {
	  vi = v[i];
	  d2 = (xg - x[i]) * (xg - x[i]) + (yg - y[i]) * (yg - y[i]);
	  w = 1.0/d2;
	  sumw  += w;
	  sumw2 += w * w;
	  sumwv += w * vi;
	  delta = vi - runmean;
	  epsilon = delta * w / sumw;
	  runmean += epsilon;
	  m2 += (sumw - w) * delta * epsilon;
	}
	/* compute ratio */
	MAT(num, ix, iy, Ny) = sumwv;
	MAT(den, ix, iy, Ny) = sumw;
	MAT(rat, ix, iy, Ny) = sumwv/sumw;
	MAT(mtwo, ix, iy, Ny) = m2;
	MAT(wtwo, ix, iy, Ny) = sumw2;
      }
    }
  } else {
    /* general case */
    for(ix = 0, xg=x0; ix < Nx; ix++, xg+=dx) {
      if(ix % 256 == 0) R_CheckUserInterrupt();
      for(iy = 0, yg=y0; iy < Ny; iy++, yg+=dy) {
	sumwv = sumw = sumw2 = m2 = runmean = 0.0;
	/* loop over data points, accumulating numerator and denominator */
	for(i = 0; i < N; i++) {
	  d2 = (xg - x[i]) * (xg - x[i]) + (yg - y[i]) * (yg - y[i]);
	  w = 1.0/pow(d2, pon2);
	  sumw  += w;
	  sumw2 += w * w;
	  sumwv += w * vi;
	  delta = vi - runmean;
	  epsilon = delta * w / sumw;
	  runmean += epsilon;
	  m2 += (sumw - w) * delta * epsilon;
	}
	/* compute ratio */
	MAT(num, ix, iy, Ny) = sumwv;
	MAT(den, ix, iy, Ny) = sumw;
	MAT(rat, ix, iy, Ny) = sumwv/sumw;
	MAT(mtwo, ix, iy, Ny) = m2;
	MAT(wtwo, ix, iy, Ny) = sumw2;
      }
    }
  }
}

/* Leave-one-out IDW at data points only */

void idwloo2(x, y, v, n, power, num, den, rat, mtwo, wtwo)
     double *x, *y, *v;           /* data points and values */
     int *n;
     double *power;                   /* exponent for IDW */
     double *num, *den, *rat, *mtwo, *wtwo;   /* output vectors - initialised 0 */
{
  int N, i, j, maxchunk;
  double xi, yi, d2, w, pon2,
    vj, sumw, sumwv, sumw2, runmean, m2, delta, epsilon;
  
  N  = *n;
  pon2 = (*power)/2.0;

  if(pon2 == 1.0) {
    /* slightly faster code when power=2 */
    OUTERCHUNKLOOP(i, N, maxchunk, 16384) {
      R_CheckUserInterrupt();
      INNERCHUNKLOOP(i, N, maxchunk, 16384) {
	xi = x[i];
	yi = y[i];
	sumwv = sumw = sumw2 = m2 = runmean = 0.0;
	if(i > 0) {
	  for(j = 0; j < i; j++) {
	    vj = v[j];
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/d2;
	    sumwv += w * vj;
	    sumw  += w;
	    sumw2 += w * w;
	    delta = vj - runmean;
	    epsilon = delta * w / sumw;
	    runmean += epsilon;
	    m2 += (sumw - w) * delta * epsilon;
	  }
	}
	if(i < N-1) {
	  for(j = i+1; j < N; j++) {
	    vj = v[j];
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/d2;
	    sumwv += w * vj;
	    sumw  += w;
	    sumw2 += w * w;
	    delta = vj - runmean;
	    epsilon = delta * w / sumw;
	    runmean += epsilon;
	    m2 += (sumw - w) * delta * epsilon;
	  }
	}
	/* compute ratio */
	num[i] = sumwv;
	den[i] = sumw;
	rat[i] = sumwv/sumw;
	mtwo[i] = m2;
	wtwo[i] = sumw2;
      }
    }
  } else {
    /* general case */
    OUTERCHUNKLOOP(i, N, maxchunk, 16384) {
      R_CheckUserInterrupt();
      INNERCHUNKLOOP(i, N, maxchunk, 16384) {
	xi = x[i];
	yi = y[i];
	sumwv = sumw = sumw2 = m2 = runmean = 0.0;
	if(i > 0) {
	  for(j = 0; j < i; j++) {
	    vj = v[j];
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/pow(d2, pon2);
	    sumwv += w * vj;
	    sumw  += w;
	    sumw2 += w * w;
	    delta = vj - runmean;
	    epsilon = delta * w / sumw;
	    runmean += epsilon;
	    m2 += (sumw - w) * delta * epsilon;
	  }
	}
	if(i < N-1) {
	  for(j = i+1; j < N; j++) {
	    vj = v[j];
	    d2 = (xi - x[j]) * (xi - x[j]) + (yi - y[j]) * (yi - y[j]);
	    w = 1.0/pow(d2, pon2);
	    sumwv += w * vj;
	    sumw  += w;
	    sumw2 += w * w;
	    delta = vj - runmean;
	    epsilon = delta * w / sumw;
	    runmean += epsilon;
	    m2 += (sumw - w) * delta * epsilon;
	  }
	}
	/* compute ratio */
	num[i] = sumwv;
	den[i] = sumw;
	rat[i] = sumwv/sumw;
	mtwo[i] = m2;
	wtwo[i] = sumw2;
      }
    }
  }
}

