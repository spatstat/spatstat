
#if (1 == 0)
/*
  knngrid.h

  Code template for C functions 
  k-nearest neighbours (k > 1) of each grid point

  THE FOLLOWING CODE ASSUMES THAT POINT PATTERN (xp, yp) IS SORTED
  IN ASCENDING ORDER OF x COORDINATE

  This code is #included multiple times in knngrid.c 
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
  Either or both DIST and WHICH may be defined.

  Copyright (C) Adrian Baddeley, Jens Oehlschlagel and Rolf Turner 2000-2013
  Licence: GPL >= 2

  $Revision: 1.6 $  $Date: 2016/02/02 01:31:50 $


*/
#endif

#undef PRINTALOT

void FNAME(nx, x0, xstep,  
	   ny, y0, ystep,   /* pixel grid dimensions */
           np, xp, yp,   /* data points */
	   kmax,
	   nnd, nnwhich, 
	   huge)
     /* inputs */
     int *nx, *ny, *np;
     double *x0, *xstep, *y0, *ystep, *huge;
     double *xp, *yp;
     int *kmax;
     /* outputs */
     double *nnd;
     int *nnwhich;
     /* some inputs + outputs are not used in all functions */
{ 
  int Nxcol, Nyrow;
  int i, j, ijpos;
  int Npoints, Nk, Nk1;
  int mleft, mright, mwhich, lastmwhich, unsorted, k, k1;
  double  X0, Y0, Xstep, Ystep;
  double d2, d2minK, xj, yi, dx, dy, dx2, hu, hu2, tmp;
  double *d2min; 
#ifdef WHICH
  int *which;
  int itmp;
#endif

  Nxcol   = *nx;
  Nyrow   = *ny;
  Npoints = *np;
  Nk      = *kmax;
  hu      = *huge;
  X0      = *x0;
  Y0      = *y0;
  Xstep   = *xstep;
  Ystep   = *ystep;

  Nk1     = Nk - 1;
  hu2      = hu * hu;

  if(Npoints == 0)
    return;

  lastmwhich = mwhich = 0;

  /* 
     create space to store the nearest neighbour distances and indices
     for the current grid point
  */

  d2min = (double *) R_alloc((size_t) Nk, sizeof(double));
#ifdef WHICH
  which = (int *) R_alloc((size_t) Nk, sizeof(int));
#endif

  /* loop over pixels */

  for(j = 0, xj = X0; j < Nxcol; j++, xj += Xstep) {

    R_CheckUserInterrupt();
    
#ifdef PRINTALOT
    Rprintf("j=%d, xj=%lf\n", j, xj); 
#endif

    for(i = 0, yi = Y0; i < Nyrow; i++, yi += Ystep) {

#ifdef PRINTALOT
      Rprintf("\ti=%d, yi = %lf\n", i, yi); 
#endif

      /* initialise nn distances and indices */
      d2minK = hu2;
      for(k = 0; k < Nk; k++) {
	d2min[k] = hu2;
#ifdef WHICH
	which[k] = -1;
#endif
      }

      if(lastmwhich < Npoints) {
	/* search forward from previous nearest neighbour  */
	for(mright = lastmwhich; mright < Npoints; ++mright)
	  {
	    dx = xp[mright] - xj;
	    dx2 = dx * dx; 
#ifdef PRINTALOT
	    Rprintf("\t\t%d\n", mright);
#endif
	    if(dx2 > d2minK) /* note that dx2 >= d2minK could break too early */
	      break;
	    dy = yp[mright] - yi;
	    d2 =  dy * dy + dx2;
	    if (d2 < d2minK) {
#ifdef PRINTALOT
	    Rprintf("\t\t\tNeighbour: d2=%lf\n", d2);
#endif
	      /* overwrite last entry in list of neighbours */
	      d2min[Nk1] = d2;
	      mwhich = mright;
#ifdef WHICH
	      which[Nk1] = mright;
#endif
	      /* bubble sort */
	      unsorted = YES;
	      for(k = Nk1; unsorted && k > 0; k--) {
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
	      /* adjust maximum distance */
	      d2minK = d2min[Nk1];
#ifdef PRINTALOT
	      Rprintf("\t\t\tUpdated d2minK=%lf\n", d2minK);
	      for(k = 0; k < Nk; k++)
		Rprintf("\t\t\t\td2min[%d]=%lf\n", k, d2min[k]);
#ifdef WHICH
	      for(k = 0; k < Nk; k++) 
		Rprintf("\t\t\t\twhich[%d]=%d\n", k, which[k]);
#endif
#endif
	    }
	  }
	/* end forward search */
      }
      if(lastmwhich > 0) {
	/* search backward from previous nearest neighbour */
	for(mleft = lastmwhich - 1; mleft >= 0; --mleft)
	  {
	    dx = xj - xp[mleft];
	    dx2 = dx * dx;
#ifdef PRINTALOT
	    Rprintf("\t\t%d\n", mleft);
#endif
	    if(dx2 > d2minK) /* note that dx2 >= d2minK could break too early */
	      break;
	    dy = yp[mleft] - yi;
	    d2 =  dy * dy + dx2;
	    if (d2 < d2minK) {
#ifdef PRINTALOT
	    Rprintf("\t\t\tNeighbour: d2=%lf\n", d2);
#endif
	      /* overwrite last entry in list of neighbours */
	      mwhich = mleft;
	      d2min[Nk1] = d2;
#ifdef WHICH
	      which[Nk1] = mleft;
#endif
	      /* bubble sort */
	      unsorted = YES;
	      for(k = Nk1; unsorted && k > 0; k--) {
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
	      /* adjust maximum distance */
	      d2minK = d2min[Nk1];
#ifdef PRINTALOT
	      Rprintf("\t\t\tUpdated d2minK=%lf\n", d2minK);
	      for(k = 0; k < Nk; k++) 
		Rprintf("\t\t\t\td2min[%d]=%lf\n", k, d2min[k]);
#ifdef WHICH
	      for(k = 0; k < Nk; k++) 
		Rprintf("\t\t\t\twhich[%d]=%d\n", k, which[k]);
#endif
#endif
	    }
	  }
	/* end backward search */
      }
      /* remember index of most recently-encountered neighbour */
      lastmwhich = mwhich;
#ifdef PRINTALOT
      Rprintf("\t\tlastmwhich=%d\n", lastmwhich);
#endif
      /* copy nn distances for grid point (i, j)
	 to output array nnd[ , i, j] 
      */
      ijpos = Nk * (i + j * Nyrow);
      for(k = 0; k < Nk; k++) {
#ifdef DIST
	nnd[ijpos + k] = sqrt(d2min[k]);
#endif
#ifdef WHICH
	nnwhich[ijpos + k] = which[k] + 1;  /* R indexing */
#endif
      }
      /* end of loop over points i */
    }
  }
}


