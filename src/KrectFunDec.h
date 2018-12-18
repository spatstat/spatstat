/*

  KrectFunDec.h

  $Revision: 1.4 $     $Date: 2018/12/18 02:43:11 $

  Function declarations for Krect

  Macros: 
      FNAME     function name

      WEIGHTED  #defined for weighted version (Kinhom etc)

  +++  Copyright (C) Adrian Baddeley 2014 ++++

*/

void FNAME(width, height,
	   nxy, x, y, 
#ifdef WEIGHTED
           w,
#endif
	   nr, rmax, trimedge, 
	   doIso, doTrans, doBord, doUnco,
	   iso, trans, bnumer, bdenom, unco)
     /* input data */
     double *width, *height;   /* window is (0, width) x (0, height) */
     int    *nxy;           /* number of (x,y) points */
     double *x, *y;         /* (x,y) coordinates */
#ifdef WEIGHTED
     double *w;             /* weights (e.g. reciprocal intensities) */
#endif
     /* algorithm parameters */
     int    *nr;            /* number of r values */
     double *rmax;          /* maximum r value */
     double *trimedge;      /* maximum edge correction weight */
     int    *doIso;         /* logical: whether to do isotropic correction */
     int    *doTrans;       /* logical: whether to do translation correction */
     int    *doBord;        /* logical: whether to do border correction */
     int    *doUnco;        /* logical: whether to do uncorrected estimator */
     /* outputs */
     /* These are vectors of length nr if required, otherwise ignored */
     double *iso;           /* isotropic-corrected estimator */
     double *trans;         /* translation-corrected estimator */
     COUNTTYPE *bnumer;        /* numerator of border-corrected estimator */
     COUNTTYPE *bdenom;        /* denominator of border-corrected estimator */
     COUNTTYPE *unco;          /* uncorrected estimator */
{
  int i, j, l, ldist, lbord, M, maxchunk, N, Nr, N1, Nr1;
  double rstep, Rmax, R2max, wide, high, trim;
  double xi, yi, bdisti, bx, by, bratio;
  double dx, dy, dx2, dij, dij2,  dratio, edgetrans, edgeiso;
  double dL, dR, dD, dU, bLU, bLD, bRU, bRD, bUL, bUR, bDL, bDR;
  double aL, aR, aD, aU, cL, cR, cU, cD, extang;
  int ncor, corner;
  COUNTTYPE *numerLowAccum, *numerHighAccum, *denomAccum;
  COUNTTYPE naccum, daccum;
  double accum;
#ifdef WEIGHTED
  double wi, wj, wij;
#endif

#ifdef WEIGHTED

#define ZERO 0.0
#define WIJ wij

#else 

#define ZERO 0
#define WIJ 1

#endif
  
  N = *nxy;

  if(N == 0) 
    return;

  Nr = *nr;
  Rmax = *rmax;

  trim = *trimedge;

  N1 = N - 1;
  Nr1 = Nr - 1;
  R2max = Rmax * Rmax;
  rstep = Rmax/Nr1;

  wide = *width;
  high = *height;

  /* Allocate and initialise scratch space - for border correction,
     but do it in all cases to keep the compiler happy */

  M = (*doBord == 1) ? Nr : 1;
  numerLowAccum  = (COUNTTYPE *) R_alloc(M, sizeof(COUNTTYPE));
  numerHighAccum = (COUNTTYPE *) R_alloc(M, sizeof(COUNTTYPE));
  denomAccum     = (COUNTTYPE *) R_alloc(M, sizeof(COUNTTYPE));
  for(l = 0; l < M; l++)
    numerLowAccum[l] = numerHighAccum[l] = denomAccum[l] = ZERO;

#include "KrectV1.h"

}

#undef ZERO
#undef WIJ
