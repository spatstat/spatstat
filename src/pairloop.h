/*

  pairloop.h

  Generic code template for loop 
  collecting contributions to point x_i
  from all points x_j such that ||x_i - x_j|| <= r

  cpp variables used:

       INITIAL_I        code executed at start of 'i' loop       
       CONTRIBUTE_IJ    code executed to compute contribution from j to i
       COMMIT_I         code executed to save total contribution to i

  C variables used:
       int i, j, n, maxchunk;
       double xi, yi, dx, dy, dx2, d2, r2max;
       double *x, *y;

  $Revision: 1.5 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#ifndef CHUNKLOOP_H
#include "chunkloop.h"
#endif

#define PAIRLOOP(INITIAL_I, CONTRIBUTE_IJ, COMMIT_I)    \
  OUTERCHUNKLOOP(i, n, maxchunk, 65536) {		\
    R_CheckUserInterrupt();				\
    INNERCHUNKLOOP(i, n, maxchunk, 65536) {		\
							\
      xi = x[i];                                        \
      yi = y[i];                                        \
                                                        \
      INITIAL_I;					\
                                                        \
      if(i > 0) {					\
	for(j=i-1; j >= 0; j--) {			\
	  dx = x[j] - xi;				\
	  dx2 = dx * dx;				\
	  if(dx2 > r2max)				\
	    break;					\
	  dy = y[j] - yi;				\
	  d2 = dx2 + dy * dy;				\
	  if(d2 <= r2max) {				\
	    CONTRIBUTE_IJ;				\
	  }						\
	}						\
      }							\
                                                        \
      if(i+1 < n) {					\
	for(j=i+1; j < n; j++) {			\
	  dx = x[j] - xi;				\
	  dx2 = dx * dx;				\
	  if(dx2 > r2max)				\
	    break;					\
	  dy = y[j] - yi;				\
	  d2 = dx2 + dy * dy;				\
	  if(d2 <= r2max) {				\
	    CONTRIBUTE_IJ;				\
	  }						\
	}						\
      }							\
      COMMIT_I;						\
    }							\
  }							
