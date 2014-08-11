/*

  nn3Ddist.c

  Nearest Neighbour Distances in 3D 

  $Revision: 1.5 $     $Date: 2012/03/27 09:32:53 $

  THE FOLLOWING FUNCTIONS ASSUME THAT z IS SORTED IN ASCENDING ORDER 

  nnd3D     Nearest neighbour distances 
  nnw3D     Nearest neighbours and their distances
  nnXw3D    Nearest neighbour from one list to another
  nnXx3D    Nearest neighbour from one list to another, with overlaps

  knnd3D    k-th nearest neighbour distances
  knnw3D    k-th nearest neighbours and their distances
*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>
#include "chunkloop.h"

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

double sqrt();

/* THE FOLLOWING CODE ASSUMES THAT z IS SORTED IN ASCENDING ORDER */

void nnd3D(n, x, y, z, nnd, huge)
/* inputs */
     int *n;
     double *x, *y, *z, *huge;
     /* output */
     double *nnd;
{ 
  int npoints, i, j, maxchunk;
  double d2, d2min, xi, yi, zi, dx, dy, dz, dz2, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {
#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

      d2min = hu2;
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
	  if(dz2 > d2min) 
	    break;
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) 
	    d2min = d2;
	}
      }
    
      /* search forward */
      if(i < npoints - 1) {
	for(j = i + 1; j < npoints; ++j) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("R");
#endif
	  dz = z[j] - zi;
	  dz2 = dz * dz;
	  if(dz2 > d2min) 
	    break;
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) 
	    d2min = d2;
	}
      }
#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      nnd[i] = sqrt(d2min);
    }
  }
}


/* nnw3D: same as nnd3D, 
   but also returns id of nearest neighbour 
*/

void nnw3D(n, x, y, z, nnd, nnwhich, huge)
/* inputs */
     int *n;
     double *x, *y, *z, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int npoints, i, j, which, maxchunk;
  double d2, d2min, xi, yi, zi, dx, dy, dz, dz2, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {
      d2min = hu2;
      which = -1;
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* search backward */
      if(i > 0){
	for(j = i - 1; j >= 0; --j) {
	  dz = z[j] - zi;
	  dz2 = dz * dz;
	  if(dz2 > d2min)
	    break;
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    which = j;
	  }
	}
      }

      /* search forward */
      if(i < npoints - 1) {
	for(j = i + 1; j < npoints; ++j) {
	  dz = z[j] - zi;
	  dz2 = dz * dz;
	  if(dz2 > d2min)
	    break;
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    which = j;
	  }
	}
      }
      nnd[i] = sqrt(d2min);
      nnwhich[i] = which;
    }
  }
}


/* 
   nnXw3D:  for TWO point patterns X and Y,
   find the nearest neighbour 
   (from each point of X to the nearest point of Y)
   returning both the distance and the identifier

   Requires both patterns to be sorted in order of increasing z coord
*/

void nnXw3D(n1, x1, y1, z1, n2, x2, y2, z2, nnd, nnwhich, huge)
/* inputs */
     int *n1, *n2;
     double *x1, *y1, *z1, *x2, *y2, *z2, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int npoints1, npoints2, i, j, jwhich, lastjwhich, maxchunk;
  double d2, d2min, x1i, y1i, z1i, dx, dy, dz, dz2, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  lastjwhich = 0;

  OUTERCHUNKLOOP(i, npoints1, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints1, maxchunk, 16384) {
      d2min = hu2;
      jwhich = -1;
      x1i = x1[i];
      y1i = y1[i];
      z1i = z1[i];

      /* search backward from previous nearest neighbour */
      if(lastjwhich > 0) {
	for(j = lastjwhich - 1; j >= 0; --j) {
	  dz = z2[j] - z1i;
	  dz2 = dz * dz;
	  if(dz2 > d2min)
	    break;
	  dx = x2[j] - x1i;
	  dy = y2[j] - y1i;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = j;
	  }
	}
      }

      /* search forward from previous nearest neighbour  */
      if(lastjwhich < npoints2) {
	for(j = lastjwhich; j < npoints2; ++j) {
	  dz = z2[j] - z1i;
	  dz2 = dz * dz;
	  if(dz2 > d2min)
	    break;
	  dx = x2[j] - x1i;
	  dy = y2[j] - y1i;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = j;
	  }
	}
      }

      nnd[i] = sqrt(d2min);
      nnwhich[i] = jwhich;
      lastjwhich = jwhich;
    }
  }
}


/* 
   nnXx3D:  similar to nnXw3D
   but allows X and Y to include common points
   (which are not to be counted as neighbours)

   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].

   Requires both patterns to be sorted in order of increasing y coord
*/

void nnXx3D(n1, x1, y1, z1, id1, n2, x2, y2, z2, id2, nnd, nnwhich, huge)
/* inputs */
     int *n1, *n2, *id1, *id2;
     double *x1, *y1, *z1, *x2, *y2, *z2, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int npoints1, npoints2, i, j, jwhich, lastjwhich, id1i, maxchunk;
  double dmin, d2, d2min, x1i, y1i, z1i, dx, dy, dz, dz2, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  lastjwhich = 0;

  for(i = 0; i < npoints1; i++) {
    
    R_CheckUserInterrupt();
    
    dmin = hu;
    d2min = hu2;
    jwhich = -1;
    x1i = x1[i];
    y1i = y1[i];
    z1i = z1[i];
    id1i = id1[i];

    /* search backward from previous nearest neighbour */
    if(lastjwhich > 0) {
      for(j = lastjwhich - 1; j >= 0; --j) {
	dz = z2[j] - z1i;
	dz2 = dz * dz;
	if(dz2 > d2min)
	  break;
	/* do not compare identical points */
	if(id2[j] != id1i) {
	  dx = x2[j] - x1i;
	  dy = y2[j] - y1i;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = j;
	  }
	}
      }
    }

    /* search forward from previous nearest neighbour  */
    if(lastjwhich < npoints2) {
      for(j = lastjwhich; j < npoints2; ++j) {
	dz = z2[j] - z1i;
	dz2 = dz * dz;
	if(dz2 > d2min)
	  break;
	/* do not compare identical points */
	if(id2[j] != id1i) {
	  dx = x2[j] - x1i;
	  dy = y2[j] - y1i;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = j;
	  }
	}
      }
    }
    nnd[i] = sqrt(d2min);
    nnwhich[i] = jwhich;
    lastjwhich = jwhich;
  }
}

/* 
   knnd3D

   nearest neighbours 1:kmax

*/

void knnd3D(n, kmax, x, y, z, nnd, huge)
/* inputs */
     int *n, *kmax;
     double *x, *y, *z, *huge;
     /* output matrix (npoints * kmax) in ROW MAJOR order */
     double *nnd;
{ 
  int npoints, nk, nk1, i, j, k, k1, unsorted, maxchunk;
  double d2, d2minK, xi, yi, zi, dx, dy, dz, dz2, hu, hu2, tmp;
  double *d2min;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the nearest neighbour distances
     for the current point
  */

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));

  /* loop over points */

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {
#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

      /* initialise nn distances */

      d2minK = hu2;
      for(k = 0; k < nk; k++) 
	d2min[k] = hu2;

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
	    /* bubble sort */
	    unsorted = TRUE;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
	      } else {
		unsorted = FALSE;
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
	    /* bubble sort */
	    unsorted = TRUE;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
	      } else {
		unsorted = FALSE;
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

      /* compute nn distances for point i 
	 and copy to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
	nnd[nk * i + k] = sqrt(d2min[k]);
      }
    }
  }
}

/* 
   knnw3D

   nearest neighbours 1:kmax

   returns distances and indices

*/

void knnw3D(n, kmax, x, y, z, nnd, nnwhich, huge)
/* inputs */
     int *n, *kmax;
     double *x, *y, *z, *huge;
     /* output matrices (npoints * kmax) in ROW MAJOR order */
     double *nnd;
     int    *nnwhich;
{ 
  int npoints, nk, nk1, i, j, k, k1, unsorted, itmp, maxchunk;
  double d2, d2minK, xi, yi, zi, dx, dy, dz, dz2, hu, hu2, tmp;
  double *d2min; 
  int *which;

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
  which = (int *) R_alloc((size_t) nk, sizeof(int));

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
	which[k] = -1;
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
	    which[nk1] = j;
	    /* bubble sort */
	    unsorted = TRUE;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
		itmp = which[k1];
		which[k1] = which[k];
		which[k] = itmp;
	      } else {
		unsorted = FALSE;
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
	    which[nk1] = j;
	    /* bubble sort */
	    unsorted = TRUE;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
		itmp = which[k1];
		which[k1] = which[k];
		which[k] = itmp;
	      } else {
		unsorted = FALSE;
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
	nnd[nk * i + k] = sqrt(d2min[k]);
	nnwhich[nk * i + k] = which[k];
      }
	
    }
  }
}

