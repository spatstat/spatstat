  /* 
     KrectBody.h 

     +++ Copyright (C) Adrian Baddeley, Julian Gilbey and Rolf Turner 2014 ++++

     Main function body for 'Krect' 

     Included multiple times with different values of the macros: 
            (#define or #undef)
     WEIGHTED
     ISOTROPIC
     TRANSLATION
     BORDER
     UNCORRECTED

     **Assumes point pattern is sorted in increasing order of x coordinate**
     **Assumes window is (0,wide) x (0, high) **
     **Assumes output vectors were initialised to zero**

     Variables are declared in 'KrectFunDec.c'

     This algorithm is optimal (amongst the choices in spatstat)
     when the window is a rectangle *and* at least one of
     the ISOTROPIC, TRANSLATION corrections is needed.
     There are faster algorithms for the border correction on its own.

     $Revision: 1.3 $ $Date: 2014/02/09 03:01:27 $

  */

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < N) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > N) maxchunk = N;

    /* ............. LOOP OVER i ................. */

    for(; i < maxchunk; i++) {

      xi = x[i];
      yi = y[i];
#ifdef WEIGHTED
      wi = w[i];
#endif

#ifdef BORDER
      /* For border correction */
      /* compute distance to border */
      bx = MIN(xi, (wide - xi));
      by = MIN(yi, (high - yi));
      bdisti = MIN(bx, by);
      /* denominator will ultimately be incremented for all r < b[i] */
      bratio = bdisti/rstep;
      /* lbord is the largest integer STRICTLY less than bratio */
      lbord = (int) ceil(bratio) - 1;
      lbord = (lbord <= Nr1) ? lbord : Nr1;
      /* increment entry corresponding to r = b[i] */
#ifdef WEIGHTED
      if(lbord >= 0) 
	denomAccum[lbord] += wi;
#else
      if(lbord >= 0) 
	(denomAccum[lbord])++;
#endif
#endif

#ifdef ISOTROPIC
      /* For isotropic correction */
      /* 
	 perpendicular distance from point i to each edge of rectangle
	 L = left, R = right, D = down, U = up
      */
      dL = xi;
      dR = wide - xi;
      dD = yi;
      dU = high - yi;
      /*
	test for corner of the rectangle
      */
      ncor = SMALL(dL) + SMALL(dR) + SMALL(dD) + SMALL(dU);
      corner = (ncor >= 2);
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
#endif

      /* ............. LOOP OVER j ................. */
      /* scan through points (x[j],y[j]) */

      /* 
	 scan backward from i-1 
	 until |x[j]-x[i]| > Rmax
      */
      if(i > 0) {
	for(j=i-1; j >= 0; j--) {
	  /* squared interpoint distance */
	  dx = xi - x[j];
	  dx2 = dx * dx;
	  if(dx2 >= R2max)
	    break;
	  dy = y[j] - yi;
	  dij2 = dx2 + dy * dy;
	  if(dij2 < R2max) {
#include "KrectIncrem.h"	    
	  }
	}
      }

      /* 
	 scan forward from i+1 
	 until x[j]-x[i] > Rmax

      */
      if(i < N1) {
	for(j=i+1; j < N; j++) {
	  /* squared interpoint distance */
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 >= R2max) 
	    break;
	  dy = y[j] - yi;
	  dij2 = dx2 + dy * dy;
	  if(dij2 < R2max) {
#include "KrectIncrem.h"	    
	  }
	}
      }
    }
  }

  /* 
    ..................  END OF LOOPS ................................
  */

  /* ............. compute cumulative functions ..................... */

#ifdef UNCORRECTED
  naccum = ZERO;
  for(l = 0; l < Nr; l++) {
    unco[l] += naccum;
    naccum = unco[l];
  }
#endif    

#ifdef ISOTROPIC
  accum = 0.0;
  for(l = 0; l < Nr; l++) {
    iso[l] += accum;
    accum = iso[l];
  }
#endif    
   
#ifdef TRANSLATION
  accum = 0.0;
  for(l = 0; l < Nr; l++) {
    trans[l] += accum;
    accum = trans[l];
  }
#endif    
   
#ifdef BORDER
  /* 
     Now use the accumulated values to compute the numerator and denominator.
     The value of denomAccum[l] should be added to denom[k] for all k <= l.
     numerHighAccum[l] should be added to numer[k] for all k <=l
     numerLowAccum[l] should then be subtracted from  numer[k] for k <= l.
  */
  for(l=Nr1, naccum=daccum=ZERO; l>=0; l--) {
    daccum += denomAccum[l];
    bdenom[l] = daccum;

    naccum += numerHighAccum[l];
    bnumer[l] = naccum;
    naccum -= numerLowAccum[l];
  }

#endif

