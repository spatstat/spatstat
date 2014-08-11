/*
  
  Knone.h

  Code template for K function estimators in Knone.c

  Variables:

  FNAME        function name

  OUTTYPE      storage type of the output 'numer' 
  ('int' or 'double')

  WEIGHTED     #defined for weighted (inhom) K function

  Copyright (C) Adrian Baddeley, Julian Gilbey and Rolf Turner 2000-2013
  Licence: GPL >= 2

  $Revision: 1.6 $     $Date: 2013/09/18 04:08:26 $

*/

void FNAME(
	   nxy, x, y, 
#ifdef WEIGHTED
	   w,
#endif
	   nr, rmax, numer) 
/* inputs */
     int *nxy, *nr;
     double *x, *y, *rmax;
#ifdef WEIGHTED
     double *w;
#endif
     /* output */
     OUTTYPE *numer;
{
  int i, j, l, n, nt, n1, lmin, lmax, maxchunk;
  double dt, tmax, tmax2, xi, yi;
  double dratio, dij, dij2, dx, dy, dx2;
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
  lmax = nt - 1;

  dt = (*rmax)/(nt-1);
  tmax = *rmax;
  tmax2 = tmax * tmax;

  /* initialise */
  for(l = 0; l < nt; l++)
    numer[l] =  ZERO;

  if(n == 0) 
    return;

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < n) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > n) maxchunk = n;

    for(; i < maxchunk; i++) {

#ifdef WEIGHTED
      wi = w[i];
#endif
      xi = x[i];
      yi = y[i];

      /* 
	 scan backward from i-1 
	 until x[j] < x[i] -tmax or until we run out 
      */
      if(i > 0) {
	for(j=i-1; j >= 0; j--) {
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 >= tmax2)
	    break;
	  dy = y[j] - yi;
	  dij2 = dx2 + dy * dy;
	  if(dij2 < tmax2) {
#ifdef WEIGHTED 
	    wj = w[j];
#endif
	    /* increment numerator for all r >= dij */
	    dij = (double) sqrt(dij2);
	    dratio = dij/dt;
	    /* smallest integer greater than or equal to dratio */
	    lmin = (int) ceil(dratio);
	    /* effectively increment entries lmin to lmax inclusive */
	    if(lmin <= lmax) {
#ifdef WEIGHTED
	      wij = wi * wj;
#endif
	      numer[lmin] += WIJ;
	    }
	  }
	}
      }

      /* 
	 scan forward from i+1 
	 until x[j] > x[i] + tmax or until we run out 
      */
      if(i < n1) {
	for(j=i+1; j < n; j++) {
	  /* squared interpoint distance */
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 >= tmax2)
	    break;
	  dy = y[j] - yi;
	  dij2 = dx2 + dy * dy;
	  if(dij2 < tmax2) {
#ifdef WEIGHTED 
	    wj = w[j];
#endif
	    /* increment numerator for all r >= dij */
	    dij = (double) sqrt(dij2);
	    dratio = dij/dt;
	    /* smallest integer greater than or equal to dratio */
	    lmin = (int) ceil(dratio);
	    /* increment entries lmin to lmax inclusive */
	    if(lmin <= lmax) {
#ifdef WEIGHTED
	      wij = wi * wj;
#endif
	      numer[lmin] += WIJ;
	    }
	  }
	}
      }
    }
  }
  /* 
     Now accumulate the numerator.
  */

  if(nt > 1)
    for(l=1; l < nt; l++)
      numer[l] += numer[l-1];

}

#undef ZERO
#undef WI 
#undef WJ 
#undef WIJ

