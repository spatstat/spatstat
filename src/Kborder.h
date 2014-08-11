/*
  
  Kborder.h

  Code template for K function estimators in Kborder.c

  Variables:

     FNAME        function name

     OUTTYPE      storage type of the output vectors
                  ('int' or 'double')

     WEIGHTED     #defined for weighted (inhom) K function


  $Revision: 1.7 $     $Date: 2012/03/18 11:29:43 $

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
  double dt, tmax, tmax2, xi, yi, bi, bi2, maxsearch, max2search;
  double bratio, dratio, dij, dij2, dx, dy, dx2;
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
  for(l = 0; l < nt; l++)
    numer[l] = denom[l] = ZERO;

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
      lmax = (int) floor(bratio);
      lup = (int) ceil(bratio);
      if(lmax == lup) --lmax;
      /* lmax is the largest integer STRICTLY less than bratio */
      lmax = (lmax <= nt1) ? lmax : nt1;
      /* increment entries 0 to lmax */
      if(lmax >= 0) {
	for(l = 0; l <= lmax; l++) 
	  denom[l] += WI;
      }

      /*  ----------  NUMERATOR -----------*/
      /* scan through points (x[j],y[j]) */
      bi2 = bi * bi;
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
	  if(dx2 > max2search)
	    break;
	  dy = y[j] - yi;
	  dij2 = dx2 + dy * dy;
	  if(dij2 < tmax2 && dij2 < bi2) {
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
	      for(l = lmin; l <= lmax; l++) 
		numer[l] += WIJ;
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
	  if(dx2 > max2search) 
	    break;
	  dy = y[j] - yi;
	  dij2 = dx2 + dy * dy;
	  if(dij2 < tmax2 && dij2 < bi2) {
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
	      for(l = lmin; l <= lmax; l++) 
		numer[l] += WIJ;
	    }
	  }
	}
      }
    }
  }
}

#undef ZERO
#undef WI 
#undef WJ
#undef WIJ

