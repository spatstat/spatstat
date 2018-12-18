/*

  knn3Ddist.h

  Code template for k-nearest-neighbour algorithms for 3D point patterns

  Input is a single point pattern - supports 'nndist' and 'nnwhich'

  This code is #included multiple times in nn3Ddist.c
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
  Either or both DIST and WHICH may be defined.

  THE FOLLOWING CODE ASSUMES THAT THE POINT PATTERN IS SORTED
  IN ASCENDING ORDER OF THE z COORDINATE

  $Revision: 1.4 $ $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

void FNAME(n, kmax, x, y, z, nnd, nnwhich, huge)
/* inputs */
     int *n, *kmax;
     double *x, *y, *z, *huge;
     /* output matrices (npoints * kmax) in ROW MAJOR order */
     double *nnd;
     int    *nnwhich;
{ 
  int npoints, nk, nk1, i, j, k, k1, unsorted, maxchunk;
  double d2, d2minK, xi, yi, zi, dx, dy, dz, dz2, hu, hu2, tmp;
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

  /* loop over points */

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {

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
      zi = z[i];

      /* search backward */
      if(i > 0) {
	for(j = i - 1; j >= 0; --j) {

#ifdef SPATSTAT_DEBUG
	  Rprintf("L");
#endif
	  dz = z[j] - zi;
	  dz2 = dz * dz; 
	  if(dz2 > d2minK)
	    break;
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2minK) {
	    /* overwrite last entry */
	    d2min[nk1] = d2;
#ifdef WHICH
	    which[nk1] = j;
#endif
	    /* bubble sort */
	    unsorted = YES;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
#ifdef WHICH
		itmp = which[k1];
		which[k1] = which[k];
		which[k] = itmp;
#endif
	      } else {
		unsorted = NO;
	      }
	    }
	    /* adjust maximum distance */
	    d2minK = d2min[nk1];
	  }
	}
      }

      /* search forward */
      if(i + 1 < npoints) {
	for(j = i + 1; j < npoints; ++j) {

#ifdef SPATSTAT_DEBUG
	  Rprintf("R");
#endif
	  dz = z[j] - zi;
	  dz2 = dz * dz;
	  if(dz2 > d2minK)
	    break;
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2minK) {
	    /* overwrite last entry */
	    d2min[nk1] = d2;
#ifdef WHICH
	    which[nk1] = j;
#endif
	    /* bubble sort */
	    unsorted = YES;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
#ifdef WHICH
		itmp = which[k1];
		which[k1] = which[k];
		which[k] = itmp;
#endif
	      } else {
		unsorted = NO;
	      }
	    }
	    /* adjust maximum distance */
	    d2minK = d2min[nk1];
	  }
	}
      }

#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      /* calculate nn distances for point i 
	 and copy to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
#ifdef DIST
	nnd[nk * i + k] = sqrt(d2min[k]);
#endif
#ifdef WHICH
	/* convert from C to R indexing */
	nnwhich[nk * i + k] = which[k] + 1;
#endif
      }
	
    }
  }
}

