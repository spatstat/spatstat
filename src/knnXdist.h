
#if (1 == 0)
/*
  knnXdist.h

  Code template for C functions supporting nncross 
  for k-nearest neighbours (k > 1)

  THE FOLLOWING CODE ASSUMES THAT LISTS ARE SORTED
  IN ASCENDING ORDER OF y COORDINATE

  This code is #included multiple times in knndistance.c 
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
        EXCLUDE   #defined if exclusion mechanism is used
  Either or both DIST and WHICH may be defined.

  When EXCLUDE is defined,
  code numbers id1, id2 are attached to the patterns X and Y respectively, 
  such that
  x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].

  Copyright (C) Adrian Baddeley, Jens Oehlschlagel and Rolf Turner 2000-2013
  Licence: GPL >= 2

  $Revision: 1.10 $  $Date: 2013/12/10 03:29:55 $


*/
#endif

void FNAME(n1, x1, y1, id1, 
           n2, x2, y2, id2, 
	   kmax,
	   nnd, nnwhich, 
	   huge)
     /* inputs */
     int *n1, *n2;
     double *x1, *y1, *x2, *y2, *huge;
     int *id1, *id2;
     int *kmax;
     /* outputs */
     double *nnd;
     int *nnwhich;
     /* some inputs + outputs are not used in all functions */
{ 
  int npoints1, npoints2, nk, nk1;
  int maxchunk, i, jleft, jright, jwhich, lastjwhich, unsorted, k, k1;
  double d2, d2minK, x1i, y1i, dx, dy, dy2, hu, hu2, tmp;
  double *d2min; 
#ifdef WHICH
  int *which;
  int itmp;
#endif
#ifdef EXCLUDE
  int id1i;
#endif
#ifdef TRACER
  int kk;
#endif

  npoints1 = *n1;
  npoints2 = *n2;
  nk       = *kmax;
  nk1      = nk - 1;
  hu       = *huge;
  hu2      = hu * hu;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  lastjwhich = 0;

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
  while(i < npoints1) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > npoints1) maxchunk = npoints1;

    for(; i < maxchunk; i++) {

      /* initialise nn distances and indices */
      d2minK = hu2;
      jwhich = -1;
      for(k = 0; k < nk; k++) {
	d2min[k] = hu2;
#ifdef WHICH
	which[k] = -1;
#endif
      }

      x1i = x1[i];
      y1i = y1[i];
#ifdef EXCLUDE
      id1i = id1[i];
#endif

#ifdef TRACER
      Rprintf("i=%d : (%lf, %lf) ..................... \n", i, x1i, y1i);
#endif

      if(lastjwhich < npoints2) {
#ifdef TRACER
	Rprintf("\tForward search from lastjwhich=%d:\n", lastjwhich);
#endif
	/* search forward from previous nearest neighbour  */
	for(jright = lastjwhich; jright < npoints2; ++jright)
	  {
#ifdef TRACER
	    Rprintf("\tjright=%d \t (%lf, %lf)\n", 
		    jright, x2[jright], y2[jright]);
#endif

	    dy = y2[jright] - y1i;
	    dy2 = dy * dy; 
#ifdef TRACER
	    Rprintf("\t\t dy2=%lf,\t d2minK=%lf\n", dy2, d2minK);
#endif
	    if(dy2 > d2minK) /* note that dy2 >= d2minK could break too early */
	      break;

#ifdef EXCLUDE
	    /* do not compare identical points */
	    if(id2[jright] != id1i) {
#ifdef TRACER
	      Rprintf("\t\t %d and %d are not identical\n", i, jright);
#endif
#endif
	      dx = x2[jright] - x1i;
	      d2 =  dx * dx + dy2;
#ifdef TRACER
	      Rprintf("\t\t d2=%lf\n", d2);
#endif
	      if (d2 < d2minK) {
		/* overwrite last entry in list of neighbours */
#ifdef TRACER
		Rprintf("\t\t overwrite d2min[nk1]=%lf by d2=%lf\n", 
			d2min[nk1], d2);
#endif
		d2min[nk1] = d2;
		jwhich = jright;
#ifdef WHICH
		which[nk1] = jright;
#endif
		/* bubble sort */
		unsorted = YES;
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
		    unsorted = NO;
		  }
		}
#ifdef TRACER
		Rprintf("\t\t sorted nn distances:\n");
		for(kk = 0; kk < nk; kk++) 
		  Rprintf("\t\t d2min[%d] = %lf\n", 
			  kk, d2min[kk]);
#endif
		/* adjust maximum distance */
		d2minK = d2min[nk1];
#ifdef TRACER
		Rprintf("\t\t d2minK=%lf\n", d2minK);
#endif
	      }
#ifdef EXCLUDE
	    }
#endif
	  }
	/* end forward search */
#ifdef TRACER
	Rprintf("\tEnd forward search\n");
#endif
      }
      if(lastjwhich > 0) {
#ifdef TRACER
	Rprintf("\tBackward search from lastjwhich=%d:\n", lastjwhich);
#endif
	/* search backward from previous nearest neighbour */
	for(jleft = lastjwhich - 1; jleft >= 0; --jleft)
	  {
#ifdef TRACER
	    Rprintf("\tjleft=%d \t (%lf, %lf)\n", 
		    jleft, x2[jleft], y2[jleft]);
#endif
	    dy = y1i - y2[jleft];
	    dy2 = dy * dy;
#ifdef TRACER
	    Rprintf("\t\t dy2=%lf,\t d2minK=%lf\n", dy2, d2minK);
#endif
	    if(dy2 > d2minK) /* note that dy2 >= d2minK could break too early */
	      break;
#ifdef EXCLUDE
	    /* do not compare identical points */
	    if(id2[jleft] != id1i) {
#ifdef TRACER
	      Rprintf("\t\t %d and %d are not identical\n", i, jleft);
#endif
#endif
	      dx = x2[jleft] - x1i;
	      d2 =  dx * dx + dy2;
#ifdef TRACER
	      Rprintf("\t\t d2=%lf\n", d2);
#endif
	      if (d2 < d2minK) {
		/* overwrite last entry in list of neighbours */
#ifdef TRACER
		Rprintf("\t\t overwrite d2min[nk1]=%lf by d2=%lf\n", 
			d2min[nk1], d2);
#endif
		d2min[nk1] = d2;
		jwhich = jleft;
#ifdef WHICH
		which[nk1] = jleft;
#endif
		/* bubble sort */
		unsorted = YES;
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
		    unsorted = NO;
		  }
		}
#ifdef TRACER
		Rprintf("\t\t sorted nn distances:\n");
		for(kk = 0; kk < nk; kk++) 
		  Rprintf("\t\t d2min[%d] = %lf\n", 
			  kk, d2min[kk]);
#endif
		/* adjust maximum distance */
		d2minK = d2min[nk1];
#ifdef TRACER
		Rprintf("\t\t d2minK=%lf\n", d2minK);
#endif
	      }
#ifdef EXCLUDE
	    }
#endif
	  }
	/* end backward search */
#ifdef TRACER
	Rprintf("\tEnd backward search\n");
#endif
      }
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
      /* save index of last neighbour encountered */
      lastjwhich = jwhich;
      /* end of loop over points i */
    }
  }
}

