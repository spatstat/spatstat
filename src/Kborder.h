/*
  
  Kborder.h

  Code template for K function estimators in Kborder.c

  Variables:

     FNAME        function name

     OUTTYPE      storage type of the output vectors
                  ('int' or 'double')

     WEIGHTED     #defined for weighted (inhom) K function


  Copyright (C) Adrian Baddeley, Julian Gilbey and Rolf Turner 2000-2013
  Licence: GPL >= 2

  $Revision: 1.9 $     $Date: 2013/04/12 06:36:00 $

*/

void FNAME(
	   nxy, x, y, 
#ifdef WEIGHTED
	   w,
#endif
	   b, nr, rmax, numer, denom) 
     /* inputs */
     int *nxy, *nr;
     double *x, *y, *b, *rmax;
#ifdef WEIGHTED
     double *w;
#endif
     /* outputs */
     OUTTYPE *numer, *denom;
{
  int i, j, l, n, nt, n1, nt1, lmin, lmax, lup, maxchunk;
  double dt, tmax, tmax2, xi, yi, bi, maxsearch, max2search;
  double bratio, dratio, dij, dij2, dx, dy, dx2;
  OUTTYPE *numerLowAccum, *numerHighAccum, *denomAccum;
  OUTTYPE naccum, daccum;
#ifdef WEIGHTED
  double wi, wj, wij;
#endif

#ifdef WEIGHTED

#define ZERO 0.0
#define WI wi
#define WJ wj
#define WIJ wij

#else 

#define ZERO 0
#define WI 1
#define WJ 1
#define WIJ 1

#endif

  n = *nxy;
  nt = *nr;

  n1 = n - 1;
  nt1 = nt - 1;

  dt = (*rmax)/(nt-1);
  tmax = *rmax;
  tmax2 = tmax * tmax;

  /* initialise */
  numerLowAccum  = (OUTTYPE *) R_alloc(nt, sizeof(OUTTYPE));
  numerHighAccum = (OUTTYPE *) R_alloc(nt, sizeof(OUTTYPE));
  denomAccum     = (OUTTYPE *) R_alloc(nt, sizeof(OUTTYPE));
  for(l = 0; l < nt; l++)
    numer[l] = denom[l] = 
      numerLowAccum[l] = numerHighAccum[l] = 
      denomAccum[l] = ZERO;

  if(n == 0) 
    return;

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < n) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > n) maxchunk = n;

    for(; i < maxchunk; i++) {

      /*  --------   DENOMINATOR  -------------*/
      bi = b[i];
#ifdef WEIGHTED
      wi = w[i];
#endif
      /* increment denominator for all r < b[i] */
      bratio = bi/dt;
      /* lmax is the largest integer STRICTLY less than bratio */
      lmax = (int) ceil(bratio) - 1;
      lmax = (lmax <= nt1) ? lmax : nt1;
      /* effectively increment entries 0 to lmax */
      if(lmax >= 0) 
	denomAccum[lmax] += WI;

      /*  ----------  NUMERATOR -----------*/
      /* scan through points (x[j],y[j]) */
      xi = x[i];
      yi = y[i];
      maxsearch = (bi < tmax) ? bi : tmax;
      max2search = maxsearch * maxsearch;

      /* 
	 scan backward from i-1 
	 until |x[j]-x[i]| > maxsearch  or until we run out 
      */
      if(i > 0) {
	for(j=i-1; j >= 0; j--) {
	  /* squared interpoint distance */
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 >= max2search)
	    break;
	  dy = y[j] - yi;
	  dij2 = dx2 + dy * dy;
	  if(dij2 < max2search) {
#ifdef WEIGHTED 
	    wj = w[j];
#endif
	    /* increment numerator for all r such that dij <= r < bi */
	    dij = (double) sqrt(dij2);
	    dratio = dij/dt;
	    /* smallest integer greater than or equal to dratio */
	    lmin = (int) ceil(dratio);
	    /* increment entries lmin to lmax inclusive */
	    if(lmax >= lmin) {
#ifdef WEIGHTED
	      wij = wi * wj;
#endif
	      numerLowAccum[lmin] += WIJ;
	      numerHighAccum[lmax] += WIJ;
	    }
	  }
	}
      }

      /* 
	 scan forward from i+1 
	 until x[j]-x[i] > maxsearch  or until we run out 

      */
      if(i < n1) {
	for(j=i+1; j < n; j++) {
	  /* squared interpoint distance */
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 >= max2search) 
	    break;
	  dy = y[j] - yi;
	  dij2 = dx2 + dy * dy;
	  if(dij2 < max2search) {
#ifdef WEIGHTED 
	    wj = w[j];
#endif
	    /* increment numerator for all r such that dij <= r < bi */
	    dij = (double) sqrt(dij2);
	    dratio = dij/dt;
	    /* smallest integer greater than or equal to dratio */
	    lmin = (int) ceil(dratio);
	    /* increment entries lmin to lmax inclusive */
	    if(lmax >= lmin) {
#ifdef WEIGHTED
	      wij = wi * wj;
#endif
	      numerLowAccum[lmin] += WIJ;
	      numerHighAccum[lmax] += WIJ;
	    }
	  }
	}
      }
    }
  }
  /* 
     Now use the accumulated values to compute the numerator and denominator.
     The value of denomAccum[l] should be added to denom[k] for all k <= l.
     numerHighAccum[l] should be added to numer[k] for all k <=l
     numerLowAccum[l] should then be subtracted from  numer[k] for k <= l.
  */

  for(l=nt1, naccum=daccum=ZERO; l>=0; l--) {
    daccum += denomAccum[l];
    denom[l] = daccum;

    naccum += numerHighAccum[l];
    numer[l] = naccum;
    naccum -= numerLowAccum[l];
  }

}

#undef ZERO
#undef WI 
#undef WJ
#undef WIJ

