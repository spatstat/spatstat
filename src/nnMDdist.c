/*

  nnMDdist.c

  Nearest Neighbour Distances in m dimensions

  $Revision: 1.18 $     $Date: 2019/10/21 11:12:32 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  Argument x is an m * n matrix 
  with columns corresponding to points
  and rows corresponding to coordinates.

  Spatial dimension m must be > 1

  THE FOLLOWING FUNCTIONS ASSUME THAT THE ROWS OF x 
  ARE SORTED IN ASCENDING ORDER OF THE FIRST COLUMN

  nndMD     Nearest neighbour distances 
  nnwMD     Nearest neighbours and their distances

  nnXwMD    Nearest neighbour from one list to another
  nnXxMD    Nearest neighbour from one list to another, with overlaps

  knndMD    k-th nearest neighbour distances
  knnwMD    k-th nearest neighbours and their distances

  knnXwMD   k-th nearest neighbours from one list to another
  knnXxMD   k-th nearest neighbours from one list to another, with overlaps

*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>
#include "chunkloop.h"

#include "yesno.h"

double sqrt();

void nndMD(n, m, x, nnd, huge)
/* inputs */
     int *n, *m;
     double *x, *huge;
     /* output */
     double *nnd;
{ 
  int npoints, mdimen, i, j, left, right, leftpos, rightpos, maxchunk;
  double d2, d2min, hu, hu2, xi0, dx0, dxj;
  double *xi;

  npoints = *n;
  mdimen  = *m; 
  xi = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

  hu = *huge;
  hu2 = hu * hu;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {

#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

      d2min = hu2;

      for(j = 0; j < mdimen; j++)
	xi[j] = x[i * mdimen + j];
      xi0 = xi[0];

#ifdef SPATSTAT_DEBUG
      Rprintf("\n (");
      for(j = 0; j < mdimen; j++)
	Rprintf("%lf, ", x[i * mdimen + j]);
      Rprintf(")\n");
#endif

    
      /* search backward */
      if(i > 0) {
	for(left = i - 1; left >= 0; --left) {

#ifdef SPATSTAT_DEBUG
	  Rprintf("L=%d, d2min=%lf\n", left, d2min);
#endif
	  dx0 = xi0 - x[left * mdimen];
	  d2 = dx0 * dx0;
	  if(d2 > d2min)
	    break;

	  leftpos = left * mdimen;
	  for(j = 1; j < mdimen && d2 < d2min; j++) {
	    dxj = xi[j] - x[leftpos + j];
	    d2 += dxj * dxj;
	  }

	  if (d2 < d2min) {
	    d2min = d2;
#ifdef SPATSTAT_DEBUG
	    Rprintf("\tupdating d2min=%lf\n", d2min);
#endif
	  }
	}
      }

      /* search forward */
      if(i < npoints - 1) {
	for(right = i + 1; right < npoints; ++right) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("R=%d, d2min=%lf\n", right, d2min);
#endif
	  dx0 = x[right * mdimen] - xi0;
	  d2  = dx0 * dx0;
	  if(d2 > d2min)
	    break;

	  rightpos = right * mdimen;
	  for(j = 1; j < mdimen && d2 < d2min; j++) {
	    dxj = xi[j] - x[rightpos + j];
	    d2 += dxj * dxj;
	  }

	  if (d2 < d2min) {
	    d2min = d2;
#ifdef SPATSTAT_DEBUG
	    Rprintf("\tupdating d2min=%lf\n", d2min);
#endif
	  }
	}
      }
#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      nnd[i] = sqrt(d2min);
    }
  }
}

/* nnwMD: same as nndMD, 
   but also returns id of nearest neighbour 
*/

void nnwMD(n, m, x, nnd, nnwhich, huge)
/* inputs */
     int *n, *m;
     double *x, *huge;
     /* output */
     double *nnd;
     int *nnwhich;
{ 
  int npoints, mdimen, i, j, left, right, leftpos, rightpos, which, maxchunk;
  double d2, d2min, hu, hu2, xi0, dx0, dxj;
  double *xi;

  npoints = *n;
  mdimen  = *m;
  xi = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

  hu = *huge;
  hu2 = hu * hu;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {

#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

      d2min = hu2;
      which = -1;

      for(j = 0; j < mdimen; j++)
	xi[j] = x[i * mdimen + j];
      xi0 = xi[0];

      /* search backward */
      if(i > 0) {
	for(left = i - 1; left >= 0; --left) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("L");
#endif
	  dx0 = xi0 - x[left * mdimen];
	  d2 = dx0 * dx0;
	  if(d2 > d2min)
	    break;
	  leftpos = left * mdimen;
	  for(j = 1; j < mdimen && d2 < d2min; j++) {
	    dxj = xi[j] - x[leftpos + j];
	    d2 += dxj * dxj;
	  }
	  if (d2 < d2min) {
	    d2min = d2;
	    which = left;
	  }
	}
      }

      /* search forward */
      if(i < npoints - 1) {
	for(right = i + 1; right < npoints; ++right) {

#ifdef SPATSTAT_DEBUG
	  Rprintf("R");
#endif
	  dx0 = x[right * mdimen] - xi0;
	  d2 = dx0 * dx0;
	  if(d2 > d2min)
	    break;

	  rightpos = right * mdimen;
	  for(j = 1; j < mdimen && d2 < d2min; j++) {
	    dxj = xi[j] - x[rightpos + j];
	    d2 += dxj * dxj;
	  }
	  if (d2 < d2min) {
	    d2min = d2;
	    which = right;
	  }
	}
      }
#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      nnd[i] = sqrt(d2min);
      /* convert index to R convention */
      nnwhich[i] = which + 1;
    }
  }
}

/* 
   nnXwMD:  for TWO point patterns X and Y,
   find the nearest neighbour 
   (from each point of X to the nearest point of Y)
   returning both the distance and the identifier

   Requires both patterns to be sorted in order of increasing first coord
*/

void nnXwMD(m, n1, x1, n2, x2, nnd, nnwhich, huge)
/* inputs */
     int *m, *n1, *n2;
     double *x1, *x2, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int mdimen, npoints1, npoints2, i, ell, jleft, jright, jwhich, lastjwhich;
  double d2, d2min, x1i0, dx0, dxell, hu, hu2;
  double *x1i;
  int maxchunk;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;
  mdimen   = *m;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  x1i = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx  = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

  lastjwhich = 0;

  OUTERCHUNKLOOP(i, npoints1, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints1, maxchunk, 16384) {
      d2min = hu2;
      jwhich = -1;
      for(ell = 0; ell < mdimen; ell++) 
	x1i[ell] = x1[i * mdimen + ell];
      x1i0 = x1i[0];

      /* search backward from previous nearest neighbour */
      if(lastjwhich > 0) {
	for(jleft = lastjwhich - 1; jleft >= 0; --jleft) {
	  dx0 = x1i0 - x2[jleft * mdimen];
	  d2 = dx0 * dx0;
	  if(d2 > d2min)
	    break;
	  for(ell = 1; ell < mdimen && d2 < d2min; ell++) {
	    dxell = x1i[ell] - x2[jleft * mdimen + ell];
	    d2 += dxell * dxell;
	  }
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = jleft;
	  }
	}
      }

      /* search forward from previous nearest neighbour  */
      if(lastjwhich < npoints2) {
	for(jright = lastjwhich; jright < npoints2; ++jright) {
	  dx0 = x2[jright * mdimen] - x1i0;
	  d2 = dx0 * dx0;
	  if(d2 > d2min) 
	    break;
	  for(ell = 1; ell < mdimen && d2 < d2min; ell++) {
	    dxell = x1i[ell] - x2[jright * mdimen + ell];
	    d2 += dxell * dxell;
	  }
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = jright;
	  }
	}
      }
      nnd[i] = sqrt(d2min);
      nnwhich[i] = jwhich + 1; /* R convention */
      lastjwhich = jwhich;
    }
  }
}


/* 
   nnXxMD:  similar to nnXwMD
   but allows X and Y to include common points
   (which are not to be counted as neighbours)

   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].

   Requires both patterns to be sorted in order of increasing first coord
*/

void nnXxMD(m, n1, x1, id1, n2, x2, id2, nnd, nnwhich, huge)
/* inputs */
     int *m, *n1, *n2;
     double *x1, *x2, *huge;
     int *id1, *id2;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int mdimen, npoints1, npoints2, i, ell, jleft, jright, jwhich, lastjwhich, id1i;
  double d2, d2min, x1i0, dx0, dxell, hu, hu2;
  double *x1i;
  int maxchunk;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;
  mdimen   = *m;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  x1i = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx  = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

  lastjwhich = 0;

  OUTERCHUNKLOOP(i, npoints1, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints1, maxchunk, 16384) {
      d2min = hu2;
      jwhich = -1;
      id1i   = id1[i];
      for(ell = 0; ell < mdimen; ell++) 
	x1i[ell] = x1[i * mdimen + ell];
      x1i0 = x1i[0];

      /* search backward from previous nearest neighbour */
      if(lastjwhich > 0) {
	for(jleft = lastjwhich - 1; jleft >= 0; --jleft) {
	  dx0 = x1i0 - x2[jleft * mdimen];
	  d2 = dx0 * dx0;
	  if(d2 > d2min)
	    break;
	  /* do not compare identical points */
	  if(id2[jleft] != id1i) {
	    for(ell = 1; ell < mdimen && d2 < d2min; ell++) {
	      dxell = x1i[ell] - x2[jleft * mdimen + ell];
	      d2 += dxell * dxell;
	    }
	    if (d2 < d2min) {
	      d2min = d2;
	      jwhich = jleft;
	    }
	  }
	}
      }

      /* search forward from previous nearest neighbour  */
      if(lastjwhich < npoints2) {
	for(jright = lastjwhich; jright < npoints2; ++jright) {
	  dx0 = x2[jright * mdimen] - x1i0;
	  d2 = dx0 * dx0;
	  if(d2 > d2min) 
	    break;
	  /* do not compare identical points */
	  if(id2[jright] != id1i) {	  
	    for(ell = 1; ell < mdimen && d2 < d2min; ell++) {
	      dxell = x1i[ell] - x2[jright * mdimen + ell];
	      d2 += dxell * dxell;
	    }
	    if (d2 < d2min) {
	      d2min = d2;
	      jwhich = jright;
	    }
	  }
	}
      }
      nnd[i] = sqrt(d2min);
      nnwhich[i] = jwhich + 1L; /* R convention */
      lastjwhich = jwhich;
    }
  }
}


/* 
   knndMD

   nearest neighbours 1:kmax

*/

void knndMD(n, m, kmax, x, nnd, huge)
/* inputs */
     int *n, *m, *kmax;
     double *x, *huge;
     /* output matrix (kmax * npoints) */
     double *nnd;
{ 
  int npoints, mdimen, nk, nk1, i, j, k, k1, left, right, unsorted, maxchunk;
  double d2, d2minK, xi0, dx0, dxj, hu, hu2, tmp;
  double *d2min, *xi;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;
  mdimen  = *m;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the squared k-th nearest neighbour distances
     for the current point
  */

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));

  /* 
     scratch space
  */
  xi = (double *) R_alloc((size_t) mdimen, sizeof(double));

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

      for(j = 0; j < mdimen; j++)
	xi[j] = x[i* mdimen + j];
      xi0 = xi[0];

#ifdef SPATSTAT_DEBUG
      Rprintf("\n (");
      for(j = 0; j < mdimen; j++)
	Rprintf("%lf, ", xi[j]);
      Rprintf(")\n");
#endif

      /* search backward */
      for(left = i - 1; left >= 0; --left) {
	dx0 = xi0 - x[left * mdimen];
	d2 = dx0 * dx0; 
	if(d2 > d2minK)
	  break;
#ifdef SPATSTAT_DEBUG
	Rprintf("L=%d\n", left);
	Rprintf("\t 0 ");
#endif
	for(j = 1; j < mdimen && d2 < d2minK; j++) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("%d ", j);
#endif
	  dxj = xi[j] - x[left * mdimen + j];
	  d2 += dxj * dxj;
	}
#ifdef SPATSTAT_DEBUG
	Rprintf("\n\t d2=%lf\n", d2);
#endif
	if (d2 < d2minK) {
	  /* overwrite last entry */
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2=%lf overwrites d2min[%d] = %lf\n", 
		  d2, nk1, d2min[nk1]);
#endif
	  d2min[nk1] = d2;
	  /* bubble sort */
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2min[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%lf, ", d2min[k]);
	  Rprintf("\n");
#endif
	  unsorted = YES;
	  for(k = nk1; unsorted && k > 0; k--) {
	    k1 = k - 1;
	    if(d2min[k] < d2min[k1]) {
	      /* swap entries */
	      tmp = d2min[k1];
	      d2min[k1] = d2min[k];
	      d2min[k] = tmp;
	    } else {
	      unsorted = NO;
	    }
	  }
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2min[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%lf, ", d2min[k]);
	  Rprintf("\n");
#endif
	  /* adjust maximum distance */
	  d2minK = d2min[nk1];
	}
      }

      /* search forward */
      for(right = i + 1; right < npoints; ++right) {

#ifdef SPATSTAT_DEBUG
	Rprintf("R=%d\n", right);
	Rprintf("\t 0 ");
#endif
	dx0 = x[right * mdimen] - xi0;
	d2 = dx0 * dx0; 
	if(d2 > d2minK)
	  break;
	for(j = 1; j < mdimen && d2 < d2minK; j++) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("%d ", j);
#endif
	  dxj = xi[j] - x[right * mdimen + j];
	  d2 += dxj * dxj;
	}
#ifdef SPATSTAT_DEBUG
	Rprintf("\n\t d2=%lf\n", d2);
#endif
	if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2=%lf overwrites d2min[%d] = %lf\n", 
		  d2, nk1, d2min[nk1]);
#endif
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  /* bubble sort */
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2min[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%lf, ", d2min[k]);
	  Rprintf("\n");
#endif
	  unsorted = YES;
	  for(k = nk1; unsorted && k > 0; k--) {
	    k1 = k - 1;
	    if(d2min[k] < d2min[k1]) {
	      /* swap entries */
	      tmp = d2min[k1];
	      d2min[k1] = d2min[k];
	      d2min[k] = tmp;
	    } else {
	      unsorted = NO;
	    }
	  }
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2min[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%lf, ", d2min[k]);
	  Rprintf("\n");
#endif
	  /* adjust maximum distance */
	  d2minK = d2min[nk1];
	}
      }

#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      /* copy nn distances for point i 
	 to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
	nnd[nk * i + k] = sqrt(d2min[k]);
      }
    }
  }
}

/* 
   knnwMD

   nearest neighbours 1:kmax

   returns distances and indices

*/


void knnwMD(n, m, kmax, x, nnd, nnwhich, huge)
/* inputs */
     int *n, *m, *kmax;
     double *x, *huge;
     /* output matrix (kmax * npoints) */
     double *nnd;
     int *nnwhich;
{ 
  int npoints, mdimen, nk, nk1, i, j, k, k1, left, right, unsorted, itmp;
  double d2, d2minK, xi0, dx0, dxj, hu, hu2, tmp;
  double *d2min, *xi;
  int *which;
  int maxchunk;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;
  mdimen  = *m;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the nearest neighbour distances and indices
     for the current point
  */

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));
  which = (int *) R_alloc((size_t) nk, sizeof(int));

  /* 
     scratch space
  */
  xi = (double *) R_alloc((size_t) mdimen, sizeof(double));

  /* loop over points */

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {

#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

      /* initialise nn distances */

      d2minK = hu2;
      for(k = 0; k < nk; k++) {
	d2min[k] = hu2;
	which[k] = -1;
      }

      for(j = 0; j < mdimen; j++)
	xi[j] = x[i* mdimen + j];
      xi0 = xi[0];

#ifdef SPATSTAT_DEBUG
      Rprintf("\n (");
      for(j = 0; j < mdimen; j++)
	Rprintf("%lf, ", x[i * mdimen + j]);
      Rprintf(")\n");
#endif

      /* search backward */
      for(left = i - 1; left >= 0; --left) {

#ifdef SPATSTAT_DEBUG
	Rprintf("L=%d, d2minK=%lf\n", left, d2minK);
	Rprintf("\t 0 ");
#endif
	dx0 = xi0 - x[left * mdimen];
	d2 = dx0 * dx0; 
	if(d2 > d2minK)
	  break;

	for(j = 1; j < mdimen && d2 < d2minK; j++) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("%d ", j);
#endif
	  dxj = xi[j] - x[left * mdimen + j];
	  d2 += dxj * dxj;
	}
	if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2=%lf overwrites d2min[%d] = %lf\n", 
		  d2, nk1, d2min[nk1]);
#endif
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  which[nk1] = left;
	  /* bubble sort */
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2min[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%lf, ", d2min[k]);
	  Rprintf("\n");
	  Rprintf("\twhich[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%d, ", which[k]);
	  Rprintf("\n");
#endif
	  unsorted = YES;
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
	      unsorted = NO;
	    }
	  }
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2min[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%lf, ", d2min[k]);
	  Rprintf("\n");
	  Rprintf("\twhich[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%d, ", which[k]);
	  Rprintf("\n");
#endif
	  /* adjust maximum distance */
	  d2minK = d2min[nk1];
	}
      }

      /* search forward */
      for(right = i + 1; right < npoints; ++right) {

#ifdef SPATSTAT_DEBUG
	Rprintf("R=%d, d2minK=%lf\n", right, d2minK);
	Rprintf("\t 0 ");
#endif
	dx0 = x[right * mdimen] - xi0;
	d2 = dx0 * dx0; 
	if(d2 > d2minK) 
	  break;
	for(j = 1; j < mdimen && d2 < d2minK; j++) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("%d ", j);
#endif
	  dxj = xi[j] - x[right * mdimen + j];
	  d2 += dxj * dxj;
	}
	if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2=%lf overwrites d2min[%d] = %lf\n", 
		  d2, nk1, d2min[nk1]);
#endif
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  which[nk1] = right;
	  /* bubble sort */
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2min[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%lf, ", d2min[k]);
	  Rprintf("\n");
	  Rprintf("\twhich[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%d, ", which[k]);
	  Rprintf("\n");
#endif
	  unsorted = YES;
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
	      unsorted = NO;
	    }
	  }
#ifdef SPATSTAT_DEBUG
	  Rprintf("\td2min[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%lf, ", d2min[k]);
	  Rprintf("\n");
	  Rprintf("\twhich[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    Rprintf("%d, ", which[k]);
	  Rprintf("\n");
#endif
	  /* adjust maximum distance */
	  d2minK = d2min[nk1];
	}
      }

#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      /* copy nn distances for point i 
	 to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
	nnd[nk * i + k] = sqrt(d2min[k]);
	/* convert index back to R convention */
	nnwhich[nk * i + k] = which[k] + 1;
      }
    }
  }
}

/* -------------- TWO POINT PATTERNS, k-nearest ------------- */

/* 
   knnXwMD

   nearest neighbours 1:kmax

   returns distances and indices

*/


void knnXwMD(m, n1, x1, n2, x2, kmax, nnd, nnwhich, huge)
  /* inputs */
  int *m, *n1, *n2, *kmax;
  double *x1, *x2, *huge;
  /* output matrix (kmax * n1) */
  double *nnd;
  int *nnwhich;
{ 
  int mdimen, npoints1, npoints2, nk, nk1;
  int i, ell, jleft, jright, jwhich, lastjwhich;
  int k, k1, unsorted, itmp;
  double d2, d2minK, x1i0, dx0, dxell, hu, hu2, tmp;
  double *d2min, *x1i;
  int *which;
  int maxchunk;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;
  mdimen   = *m;
  nk       = *kmax;
  nk1      = nk - 1;

  /* 
     create space to store the nearest neighbour distances and indices
     for the current point
  */

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));
  which = (int *) R_alloc((size_t) nk, sizeof(int));

  /* 
     scratch space for current 'from' point coordinates
  */
  x1i = (double *) R_alloc((size_t) mdimen, sizeof(double));

  lastjwhich = 0;
  
  /* loop over 'from' points */

  OUTERCHUNKLOOP(i, npoints1, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints1, maxchunk, 16384) {

#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

      /* initialise nn distances */
      d2minK = hu2;
      jwhich = -1;
      for(k = 0; k < nk; k++) {
	d2min[k] = hu2;
	which[k] = -1;
      }
      /* copy coordinates of current 'from' point */
      for(ell = 0; ell < mdimen; ell++)
	x1i[ell] = x1[i* mdimen + ell];
      x1i0 = x1i[0];

#ifdef SPATSTAT_DEBUG
      Rprintf("\n From (");
      for(ell = 0; ell < mdimen; ell++)
	Rprintf("%lf, ", x1[i * mdimen + ell]);
      Rprintf(")\n");
#endif

      if(lastjwhich > 0) {
	/* search backward from previous nearest neighbour */
	for(jleft = lastjwhich - 1; jleft >= 0; --jleft) {

#ifdef SPATSTAT_DEBUG
	  Rprintf("L=%d, d2minK=%lf\n", jleft, d2minK);
	  Rprintf("\t 0 ");
#endif
	  dx0 = x1i0 - x2[jleft * mdimen];
	  d2 = dx0 * dx0; 
	  if(d2 > d2minK)
	    break;

	  for(ell = 1; ell < mdimen && d2 < d2minK; ell++) {
#ifdef SPATSTAT_DEBUG
	    Rprintf("%d ", ell);
#endif
	    dxell = x1i[ell] - x2[jleft * mdimen + ell];
	    d2 += dxell * dxell;
	  }
	  if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	    Rprintf("\td2=%lf overwrites d2min[%d] = %lf\n", 
		    d2, nk1, d2min[nk1]);
#endif
	    /* overwrite last entry */
	    d2min[nk1] = d2;
	    which[nk1] = jleft;
	    jwhich = jleft;
	    /* bubble sort */
#ifdef SPATSTAT_DEBUG
	    Rprintf("\td2min[] before bubble sort:");
	    for(k = 0; k < nk; k++)
	      Rprintf("%lf, ", d2min[k]);
	    Rprintf("\n");
	    Rprintf("\twhich[] before bubble sort:");
	    for(k = 0; k < nk; k++)
	      Rprintf("%d, ", which[k]);
	    Rprintf("\n");
#endif
	    unsorted = YES;
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
		unsorted = NO;
	      }
	    }
#ifdef SPATSTAT_DEBUG
	    Rprintf("\td2min[] after bubble sort:");
	    for(k = 0; k < nk; k++)
	      Rprintf("%lf, ", d2min[k]);
	    Rprintf("\n");
	    Rprintf("\twhich[] after bubble sort:");
	    for(k = 0; k < nk; k++)
	      Rprintf("%d, ", which[k]);
	    Rprintf("\n");
#endif
	    /* adjust maximum distance */
	    d2minK = d2min[nk1];
	  }
	}
      }

      
      /* search forward */
      if(lastjwhich < npoints2) {
	for(jright = lastjwhich; jright < npoints2; ++jright) {

#ifdef SPATSTAT_DEBUG
	  Rprintf("R=%d, d2minK=%lf\n", jright, d2minK);
	  Rprintf("\t 0 ");
#endif
	  dx0 = x2[jright * mdimen] - x1i0;
	  d2 = dx0 * dx0; 
	  if(d2 > d2minK) 
	    break;
	  for(ell = 1; ell < mdimen && d2 < d2minK; ell++) {
#ifdef SPATSTAT_DEBUG
	    Rprintf("%d ", ell);
#endif
	    dxell = x1i[ell] - x2[jright * mdimen + ell];
	    d2 += dxell * dxell;
	  }
	  if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	    Rprintf("\td2=%lf overwrites d2min[%d] = %lf\n", 
		    d2, nk1, d2min[nk1]);
#endif
	    /* overwrite last entry */
	    d2min[nk1] = d2;
	    which[nk1] = jright;
	    jwhich = jright;
	    /* bubble sort */
#ifdef SPATSTAT_DEBUG
	    Rprintf("\td2min[] before bubble sort:");
	    for(k = 0; k < nk; k++)
	      Rprintf("%lf, ", d2min[k]);
	    Rprintf("\n");
	    Rprintf("\twhich[] before bubble sort:");
	    for(k = 0; k < nk; k++)
	      Rprintf("%d, ", which[k]);
	    Rprintf("\n");
#endif
	    unsorted = YES;
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
		unsorted = NO;
	      }
	    }
#ifdef SPATSTAT_DEBUG
	    Rprintf("\td2min[] after bubble sort:");
	    for(k = 0; k < nk; k++)
	      Rprintf("%lf, ", d2min[k]);
	    Rprintf("\n");
	    Rprintf("\twhich[] after bubble sort:");
	    for(k = 0; k < nk; k++)
	      Rprintf("%d, ", which[k]);
	    Rprintf("\n");
#endif
	    /* adjust maximum distance */
	    d2minK = d2min[nk1];
	  }
	}
      }

#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      /* copy nn distances for point i 
	 to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
	nnd[nk * i + k] = sqrt(d2min[k]);
	/* convert index back to R convention */
	nnwhich[nk * i + k] = which[k] + 1;
      }
      
      /* save index of last neighbour encountered */
      lastjwhich = jwhich;
    }
  }
}


/* 
   knnXxMD

   nearest neighbours 1:kmax with exclusions

   returns distances and indices

*/


void knnXxMD(m, n1, x1, id1, n2, x2, id2, kmax, nnd, nnwhich, huge)
  /* inputs */
  int *m, *n1, *n2, *kmax;
  double *x1, *x2, *huge;
  int *id1, *id2;
  /* output matrix (kmax * n1) */
  double *nnd;
  int *nnwhich;
{ 
  int mdimen, npoints1, npoints2, nk, nk1;
  int i, ell, jleft, jright, jwhich, lastjwhich;
  int k, k1, unsorted, itmp, id1i;
  double d2, d2minK, x1i0, dx0, dxell, hu, hu2, tmp;
  double *d2min, *x1i;
  int *which;
  int maxchunk;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;
  mdimen   = *m;
  nk       = *kmax;
  nk1      = nk - 1;

  /* 
     create space to store the nearest neighbour distances and indices
     for the current point
  */

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));
  which = (int *) R_alloc((size_t) nk, sizeof(int));

  /* 
     scratch space for current 'from' point coordinates
  */
  x1i = (double *) R_alloc((size_t) mdimen, sizeof(double));

  lastjwhich = 0;
  
  /* loop over 'from' points */

  OUTERCHUNKLOOP(i, npoints1, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints1, maxchunk, 16384) {

#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

      /* initialise nn distances */
      d2minK = hu2;
      jwhich = -1;
      for(k = 0; k < nk; k++) {
	d2min[k] = hu2;
	which[k] = -1;
      }
      /* copy coordinates of current 'from' point */
      for(ell = 0; ell < mdimen; ell++)
	x1i[ell] = x1[i* mdimen + ell];
      x1i0 = x1i[0];

      id1i = id1[i];
      
#ifdef SPATSTAT_DEBUG
      Rprintf("\n From (");
      for(ell = 0; ell < mdimen; ell++)
	Rprintf("%lf, ", x1[i * mdimen + ell]);
      Rprintf(")\n");
#endif

      if(lastjwhich > 0) {
	/* search backward from previous nearest neighbour */
	for(jleft = lastjwhich - 1; jleft >= 0; --jleft) {

#ifdef SPATSTAT_DEBUG
	  Rprintf("L=%d, d2minK=%lf\n", jleft, d2minK);
	  Rprintf("\t 0 ");
#endif
	  dx0 = x1i0 - x2[jleft * mdimen];
	  d2 = dx0 * dx0; 
	  if(d2 > d2minK)
	    break;

	  /* don't compare identical points */
	  if(id2[jleft] != id1i) {
	    for(ell = 1; ell < mdimen && d2 < d2minK; ell++) {
#ifdef SPATSTAT_DEBUG
	      Rprintf("%d ", ell);
#endif
	      dxell = x1i[ell] - x2[jleft * mdimen + ell];
	      d2 += dxell * dxell;
	    }
	    if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	      Rprintf("\td2=%lf overwrites d2min[%d] = %lf\n", 
		      d2, nk1, d2min[nk1]);
#endif
	      /* overwrite last entry */
	      d2min[nk1] = d2;
	      which[nk1] = jleft;
	      jwhich = jleft;
	      /* bubble sort */
#ifdef SPATSTAT_DEBUG
	      Rprintf("\td2min[] before bubble sort:");
	      for(k = 0; k < nk; k++)
		Rprintf("%lf, ", d2min[k]);
	      Rprintf("\n");
	      Rprintf("\twhich[] before bubble sort:");
	      for(k = 0; k < nk; k++)
		Rprintf("%d, ", which[k]);
	      Rprintf("\n");
#endif
	      unsorted = YES;
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
		  unsorted = NO;
		}
	      }
#ifdef SPATSTAT_DEBUG
	      Rprintf("\td2min[] after bubble sort:");
	      for(k = 0; k < nk; k++)
		Rprintf("%lf, ", d2min[k]);
	      Rprintf("\n");
	      Rprintf("\twhich[] after bubble sort:");
	      for(k = 0; k < nk; k++)
		Rprintf("%d, ", which[k]);
	      Rprintf("\n");
#endif
	      /* adjust maximum distance */
	      d2minK = d2min[nk1];
	    }
	  }
	}
      }

      
      /* search forward */
      if(lastjwhich < npoints2) {
	for(jright = lastjwhich; jright < npoints2; ++jright) {

#ifdef SPATSTAT_DEBUG
	  Rprintf("R=%d, d2minK=%lf\n", jright, d2minK);
	  Rprintf("\t 0 ");
#endif
	  dx0 = x2[jright * mdimen] - x1i0;
	  d2 = dx0 * dx0; 
	  if(d2 > d2minK) 
	    break;

	  /* don't compare identical points */
	  if(id2[jright] != id1i) {
	    for(ell = 1; ell < mdimen && d2 < d2minK; ell++) {
#ifdef SPATSTAT_DEBUG
	      Rprintf("%d ", ell);
#endif
	      dxell = x1i[ell] - x2[jright * mdimen + ell];
	      d2 += dxell * dxell;
	    }
	    if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	      Rprintf("\td2=%lf overwrites d2min[%d] = %lf\n", 
		      d2, nk1, d2min[nk1]);
#endif
	      /* overwrite last entry */
	      d2min[nk1] = d2;
	      which[nk1] = jright;
	      jwhich = jright;
	      /* bubble sort */
#ifdef SPATSTAT_DEBUG
	      Rprintf("\td2min[] before bubble sort:");
	      for(k = 0; k < nk; k++)
		Rprintf("%lf, ", d2min[k]);
	      Rprintf("\n");
	      Rprintf("\twhich[] before bubble sort:");
	      for(k = 0; k < nk; k++)
		Rprintf("%d, ", which[k]);
	      Rprintf("\n");
#endif
	      unsorted = YES;
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
		  unsorted = NO;
		}
	      }
#ifdef SPATSTAT_DEBUG
	      Rprintf("\td2min[] after bubble sort:");
	      for(k = 0; k < nk; k++)
		Rprintf("%lf, ", d2min[k]);
	      Rprintf("\n");
	      Rprintf("\twhich[] after bubble sort:");
	      for(k = 0; k < nk; k++)
		Rprintf("%d, ", which[k]);
	      Rprintf("\n");
#endif
	      /* adjust maximum distance */
	      d2minK = d2min[nk1];
	    }
	  }
	}
      }

#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      /* copy nn distances for point i 
	 to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
	nnd[nk * i + k] = sqrt(d2min[k]);
	/* convert index back to R convention */
	nnwhich[nk * i + k] = which[k] + 1;
      }
      
      /* save index of last neighbour encountered */
      lastjwhich = jwhich;
    }
  }
}
  
