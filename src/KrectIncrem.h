/*
  KrectIncrem.h

  Code to increment numerators of K-function

  $Revision: 1.5 $  $Date: 2014/02/09 03:00:51 $

  +++  Copyright (C) Adrian Baddeley, Julian Gilbey and Rolf Turner 2014 ++++

*/

#ifdef WEIGHTED 
	      wj = w[j];
	      wij = wi * wj;
#endif
	      /* determine index of entry to be incremented */
	      dij = (double) sqrt(dij2);
	      dratio = dij/rstep;
	      /* smallest integer greater than or equal to dratio */
	      ldist = (int) ceil(dratio);

#ifdef UNCORRECTED
	      /* ............  uncorrected estimate ................. */
#ifdef WEIGHTED
              unco[ldist] += wij;             
#else
              (unco[ldist])++;
#endif
#endif

#ifdef BORDER
	      /* ............  border correction ................. */
	      /* increment numerator for all r such that dij <= r < bi */
	      /* increment entries ldist to lbord inclusive */
#ifdef WEIGHTED
	      if(lbord >= ldist) {
		numerLowAccum[ldist] += wij;
		numerHighAccum[lbord] += wij;
	      }
#else
	      if(lbord >= ldist) {
		(numerLowAccum[ldist])++;
		(numerHighAccum[lbord])++;
	      }
#endif
#endif

#ifdef TRANSLATION
	      /* ............  translation correction ................. */
              edgetrans = 1.0/((1.0 - ABS(dx)/wide) * (1.0 - ABS(dy)/high));
              edgetrans = MIN(edgetrans, trim);
#ifdef WEIGHTED
	      trans[ldist] += wij * edgetrans;
#else
	      trans[ldist] += edgetrans;
#endif
#endif

#ifdef ISOTROPIC
	      /* ............  isotropic correction ................. */
	      /*
		half the angle subtended by the intersection between
		the circle of radius d[i,j] centred on point i
		and each edge of the rectangle (prolonged to an infinite line)
	      */
	      aL = (dL < dij) ? acos(dL/dij) : 0.0;
	      aR = (dR < dij) ? acos(dR/dij) : 0.0;
	      aD = (dD < dij) ? acos(dD/dij) : 0.0;
	      aU = (dU < dij) ? acos(dU/dij) : 0.0;

	      /* apply maxima */

	      cL = MIN(aL, bLU) + MIN(aL, bLD);
	      cR = MIN(aR, bRU) + MIN(aR, bRD);
	      cU = MIN(aU, bUL) + MIN(aU, bUR);
	      cD = MIN(aD, bDL) + MIN(aD, bDR);

	      /* total exterior angle over 2 pi */
	      extang = (cL + cR + cU + cD)/TWOPI;

	      /* add pi/2 for corners */
	      if(corner) 
		extang += 1/4;

	      /* edge correction factor */
	      edgeiso = 1 / (1 - extang);
              edgeiso = MIN(edgeiso, trim);

#ifdef WEIGHTED
	      iso[ldist] += wij * edgeiso;
#else
	      iso[ldist] += edgeiso;
#endif
#endif
