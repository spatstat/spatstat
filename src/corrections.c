/*

  corrections.c

  Edge corrections

  $Revision: 1.15 $     $Date: 2019/03/30 01:15:48 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"
#include "yesno.h"
#include "constants.h"

#undef DEBUG


/* This constant is defined in Rmath.h */
#define TWOPI M_2PI

#define MIN(A,B) (((A) < (B)) ? (A) : (B))

#define BETWEEN(X,X0,X1) ((double) ( ( (X) - (X0) ) * ( (X) - (X1) ) ) <= 0.0)

#define UNDER(X,Y,X0,Y0,X1,Y1) \
  ((double) ( ( (Y1) - (Y0) ) * ( (X) - (X0) ) ) >= (double) ( ( (Y) - (Y0) ) * ( (X1) - (X0) ) ) )

#define UNDERNEATH(X,Y,X0,Y0,X1,Y1) \
  ((((double) (X0)) < ((double) (X1))) ? UNDER(X,Y,X0,Y0,X1,Y1) : UNDER(X,Y,X1,Y1,X0,Y0))

#define TESTINSIDE(X,Y,X0,Y0,X1,Y1) \
  (BETWEEN(X,X0,X1) && UNDERNEATH(X, Y, X0, Y0, X1, Y1))


void ripleybox(nx, x, y, rmat, nr, xmin, ymin, xmax, ymax,  epsilon, out)
     /* inputs */
     int *nx, *nr;  /* dimensions */
     double *x, *y; /* coordinate vectors of length nx */
     double *rmat;  /* matrix nx by nr  */
     double *xmin, *ymin, *xmax, *ymax;  /* box dimensions */
     double *epsilon; /* threshold for proximity to corner */
     /* output */
     double *out;  /* output matrix nx by nr */
{
  int i, j, n, m, ijpos, ncor, maxchunk;
  double xx, yy, x0, y0, x1, y1, dL, dR, dU, dD, aL, aU, aD, aR, rij;
  double cL, cU, cD, cR, bLU, bLD, bRU, bRD, bUL, bUR, bDL, bDR;
  double corner, extang;
  double eps;

  n  = *nx;
  m  = *nr;
  x0 = *xmin;
  y0 = *ymin;
  x1 = *xmax;
  y1 = *ymax;
  eps = *epsilon;

  OUTERCHUNKLOOP(i, n, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, n, maxchunk, 16384) {
      xx = x[i];
      yy = y[i];
      /* 
	 perpendicular distance from point to each edge of rectangle
	 L = left, R = right, D = down, U = up
      */
      dL = xx - x0;
      dR = x1 - xx;
      dD = yy - y0;
      dU = y1 - yy;

      /*
	test for corner of the rectangle
      */
#define ABS(X) (((X) >= 0) ? (X) : (-X))
#define SMALL(X) ((ABS(X) < eps) ? 1 : 0)

      ncor = SMALL(dL) + SMALL(dR) + SMALL(dD) + SMALL(dU);
      corner = (ncor >= 2) ? YES : NO;
  
      /* 
	 angle between 
	 - perpendicular to edge of rectangle
	 and 
	 - line from point to corner of rectangle

      */
      bLU = atan2(dU, dL);
      bLD = atan2(dD, dL);
      bRU = atan2(dU, dR);
      bRD = atan2(dD, dR);
      bUL = atan2(dL, dU);
      bUR = atan2(dR, dU);
      bDL = atan2(dL, dD);
      bDR = atan2(dR, dD);

      for(j = 0; j < m; j++) {
	ijpos = j * n + i;
	rij = rmat[ijpos];
#ifdef DEBUG
	Rprintf("rij = %lf\n", rij);
#endif
	/*
	  half the angle subtended by the intersection between
	  the circle of radius r[i,j] centred on point i
	  and each edge of the rectangle (prolonged to an infinite line)
	*/
	aL = (dL < rij) ? acos(dL/rij) : 0.0;
	aR = (dR < rij) ? acos(dR/rij) : 0.0;
	aD = (dD < rij) ? acos(dD/rij) : 0.0;
	aU = (dU < rij) ? acos(dU/rij) : 0.0;
#ifdef DEBUG
	Rprintf("aL = %lf\n", aL);
	Rprintf("aR = %lf\n", aR);
	Rprintf("aD = %lf\n", aD);
	Rprintf("aU = %lf\n", aU);
#endif
	/* apply maxima */

	cL = MIN(aL, bLU) + MIN(aL, bLD);
	cR = MIN(aR, bRU) + MIN(aR, bRD);
	cU = MIN(aU, bUL) + MIN(aU, bUR);
	cD = MIN(aD, bDL) + MIN(aD, bDR);
#ifdef DEBUG
	Rprintf("cL = %lf\n", cL);
	Rprintf("cR = %lf\n", cR);
	Rprintf("cD = %lf\n", cD);
	Rprintf("cU = %lf\n", cU);
#endif

	/* total exterior angle over 2 pi */
	extang = (cL + cR + cU + cD)/TWOPI;

	/* add pi/2 for corners */
	if(corner) 
	  extang += 1/4;

#ifdef DEBUG
	Rprintf("extang = %lf\n", extang);
#endif
	/* OK, now compute weight */
	out[ijpos] = 1 / (1 - extang);
      }
    }
  }
}


/* C function ripleypoly */
#undef DEBUGPOLY
#define RIPLEYFUN ripleypoly
#include "ripleypoly.h"
#undef RIPLEYFUN

/* C function rippolDebug */
#define RIPLEYFUN rippolDebug
#define DEBUGPOLY
#include "ripleypoly.h"
#undef RIPLEYFUN
#undef DEBUGPOLY


