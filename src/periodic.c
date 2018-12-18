/* 
   periodic.c

   Routines for periodic edge correction 

   Naive algorithms O(n^2) in time (but memory-efficient)
   which can easily be adapted to more general metrics.

   Coordinates are NOT assumed to be sorted
   
   $Revision: 1.4 $ $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

#define OK 0
#define ERR_OVERFLOW 1
#define ERR_ALLOC 2

#define intRealloc(PTR, OLDLENGTH, NEWLENGTH) \
  (int *) S_realloc((char *) PTR, NEWLENGTH, OLDLENGTH, sizeof(int))

#define dblRealloc(PTR, OLDLENGTH, NEWLENGTH) \
  (double *) S_realloc((char *) PTR, NEWLENGTH, OLDLENGTH, sizeof(double))

double sqrt();

/* counterpart of 'closepairs' */

SEXP closePpair(SEXP xx,    /* spatial coordinates */
		SEXP yy,
		SEXP pp,    /* period */
		SEXP rr,    /* max distance */
		SEXP nguess) 
{
  double *x, *y;
  double xi, yi, rmax, r2max, dx, dy, d2, dxp, dyp;
  int n, k, kmax, kmaxold, maxchunk, i, j, m;
  double *period;
  double xperiod, yperiod;
  /* local storage */
  int *iout, *jout;
  double *dout;
  /* R objects in return value */
  SEXP Out, iOut, jOut, dOut;
  /* external storage pointers */
  int *iOutP, *jOutP;
  double *dOutP;

  /* protect R objects from garbage collector */
  PROTECT(xx     = AS_NUMERIC(xx));
  PROTECT(yy     = AS_NUMERIC(yy));
  PROTECT(pp     = AS_NUMERIC(pp));
  PROTECT(rr     = AS_NUMERIC(rr));
  PROTECT(nguess = AS_INTEGER(nguess));
  /* that's 5 protected arguments */
#define NINPUTS 5
  
  /* Translate arguments from R to C */
  x = NUMERIC_POINTER(xx);
  y = NUMERIC_POINTER(yy);
  n = LENGTH(xx);
  period = NUMERIC_POINTER(pp);
  xperiod = period[0];
  yperiod = period[1];
  rmax = *(NUMERIC_POINTER(rr));
  r2max = rmax * rmax;
  kmax = *(INTEGER_POINTER(nguess));

  k = 0;   /* k is the next available storage location 
              and also the current length of the list */ 

  if(n > 0 && kmax > 0) {
    /* allocate space */
    iout = (int *) R_alloc(kmax, sizeof(int));
    jout = (int *) R_alloc(kmax, sizeof(int));
    dout  =  (double *) R_alloc(kmax, sizeof(double));

    /* loop in chunks of 2^16 */
    i = 0; maxchunk = 0; 
    while(i < n) {

      R_CheckUserInterrupt();

      maxchunk += 65536; 
      if(maxchunk > n) maxchunk = n;

      for(; i < maxchunk; i++) {

	xi = x[i];
	yi = y[i];

	if(i > 0) {
	  /* scan backward */
	  for(j = i - 1; j >= 0; j--) {

	    dx = x[j] - xi;
	    if(dx < 0.0) dx = -dx;
	    dxp = xperiod - dx;
	    if(dxp < dx) dx = dxp;

	    if(dx < rmax) {

	      dy = y[j] - yi;
	      if(dy < 0.0) dy = -dy;
	      dyp = yperiod - dy;
	      if(dyp < dy) dy = dyp;

	      d2 = dx * dx + dy * dy;

	      if(d2 <= r2max) {
		/* add this (i, j) pair to output */
		if(k >= kmax) {
		  /* overflow; allocate more space */
		  kmaxold = kmax;
		  kmax    = 2 * kmax;
		  iout  = intRealloc(iout,  kmaxold, kmax);
		  jout  = intRealloc(jout,  kmaxold, kmax);
		  dout  = dblRealloc(dout,  kmaxold, kmax); 
		}
		jout[k] = j + 1; /* R indexing */
		iout[k] = i + 1;
		dout[k] = sqrt(d2);
		++k;
	      }
	    }
	  }
	}
	  
	if(i + 1 < n) {
	  /* scan forward */
	  for(j = i + 1; j < n; j++) {

	    dx = x[j] - xi;
	    if(dx < 0.0) dx = -dx;
	    dxp = xperiod - dx;
	    if(dxp < dx) dx = dxp;

	    if(dx < rmax) {

	      dy = y[j] - yi;
	      if(dy < 0.0) dy = -dy;
	      dyp = yperiod - dy;
	      if(dyp < dy) dy = dyp;

	      d2 = dx * dx + dy * dy;

	      if(d2 <= r2max) {
		/* add this (i, j) pair to output */
		if(k >= kmax) {
		  /* overflow; allocate more space */
		  kmaxold = kmax;
		  kmax    = 2 * kmax;
		  iout  = intRealloc(iout,  kmaxold, kmax);
		  jout  = intRealloc(jout,  kmaxold, kmax);
		  dout  = dblRealloc(dout,  kmaxold, kmax); 
		}
		jout[k] = j + 1; /* R indexing */
		iout[k] = i + 1;
		dout[k] = sqrt(d2);
		++k;
	      }
	    }
	  }
	}
	/* end of i loop */
      }
    }    
  }

  /* return a list of vectors */
  PROTECT(Out   = NEW_LIST(3));
  PROTECT(iOut  = NEW_INTEGER(k));
  PROTECT(jOut  = NEW_INTEGER(k));
  PROTECT(dOut  = NEW_NUMERIC(k));
#define NALLOCATED 4

  /* copy results into return object */
  if(k > 0) {
    iOutP  = INTEGER_POINTER(iOut);
    jOutP  = INTEGER_POINTER(jOut);
    dOutP  = NUMERIC_POINTER(dOut);
    for(m = 0; m < k; m++) {
      iOutP[m] = iout[m];
      jOutP[m] = jout[m];
      dOutP[m]  = dout[m];
    }
  }
  SET_VECTOR_ELT(Out, 0,  iOut);
  SET_VECTOR_ELT(Out, 1,  jOut);
  SET_VECTOR_ELT(Out, 2,  dOut);

  /* relinquish and return */
  UNPROTECT(NINPUTS+NALLOCATED); 
  return(Out);
}

