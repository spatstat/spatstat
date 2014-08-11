/*
  knndist.h

  Code template for C functions supporting knndist and knnwhich 

  THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER

  This code is #included multiple times in nndistance.c 
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
  Either or both DIST and WHICH may be defined.

  Copyright (C) Adrian Baddeley, Jens Oehlschlagel and Rolf Turner 2000-2012
  Licence: GPL >= 2

  $Revision: 1.1 $  $Date: 2012/03/18 07:08:29 $

*/

void FNAME(n, kmax, x, y, 
#ifdef DIST 
	   nnd, 
#endif
#ifdef WHICH
	   nnwhich, 
#endif
	   huge)
     /* inputs */
     int *n, *kmax;
     double *x, *y, *huge;
     /* output matrices (npoints * kmax) in ROW MAJOR order */
#ifdef DIST
     double *nnd;
#endif
#ifdef WHICH
     int    *nnwhich;
#endif
{ 
  int npoints, maxchunk, nk, nk1, i, k, k1, left, right, unsorted;
  double d2, d2minK, xi, yi, dx, dy, dy2, hu, hu2, tmp;
  double *d2min; 
#ifdef WHICH
  int *which;
  int itmp;
#endif

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the nearest neighbour distances and indices
     for the current point
  */

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));
#ifdef WHICH
  which = (int *) R_alloc((size_t) nk, sizeof(int));
#endif

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < npoints) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > npoints) maxchunk = npoints;

    for(; i < maxchunk; i++) {

#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

    /* initialise nn distances and indices */

      d2minK = hu2;
      for(k = 0; k < nk; k++) {
	d2min[k] = hu2;
#ifdef WHICH
	which[k] = -1;
#endif
      }

      xi = x[i];
      yi = y[i];

      /* search backward */
      for(left = i - 1; left >= 0; --left)
      {

#ifdef SPATSTAT_DEBUG
	Rprintf("L");
#endif
	dy = yi - y[left];
	dy2 = dy * dy;
	if(dy2 > d2minK)
	  break;

	dx = x[left] - xi;
	d2 =  dx * dx + dy2;
	if (d2 < d2minK) {
	  /* overwrite last entry */
	  d2min[nk1] = d2;
#ifdef WHICH
	  which[nk1] = left;
#endif
	  /* bubble sort */
	  unsorted = TRUE;
	  for(k = nk1; unsorted && k > 0; k--) {
	    k1 = k - 1;
	    if(d2min[k] < d2min[k1]) {
	      /* swap entries */
	      tmp  = d2min[k1];
	      d2min[k1] = d2min[k];
	      d2min[k] = tmp;
#ifdef WHICH
	      itmp = which[k1];
	      which[k1] = which[k];
	      which[k] = itmp;
#endif
	    } else {
	      unsorted = FALSE;
	    }
	  }
	  /* adjust maximum distance */
	  d2minK = d2min[nk1];
	}
      }

      /* search forward */
      for(right = i + 1; right < npoints; ++right)
	{

#ifdef SPATSTAT_DEBUG
	  Rprintf("R");
#endif
	  dy = y[right] - yi;
	  dy2 = dy * dy;
	  if(dy2 > d2minK)
	    break;

	  dx = x[right] - xi;
	  d2 =  dx * dx + dy2;
	  if (d2 < d2minK) {
	    /* overwrite last entry */
	    d2min[nk1] = d2;
#ifdef WHICH
	    which[nk1] = right;
#endif
	    /* bubble sort */
	    unsorted = TRUE;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp  = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
#ifdef WHICH
		itmp = which[k1];
		which[k1] = which[k];
		which[k] = itmp;
#endif
	      } else {
		unsorted = FALSE;
	      }
	    }
	    /* adjust maximum distance */
	    d2minK = d2min[nk1];
	  }
	}

      /* search finished for point i */

#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      /* copy nn distances for point i 
	 to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
#ifdef DIST
	nnd[nk * i + k] = sqrt(d2min[k]);
#endif
#ifdef WHICH
	nnwhich[nk * i + k] = which[k] + 1;  /* R indexing */
#endif
      }

      /* end of i loop */
    }
  }
}


