
#if (1 == 0)
/*
  nngrid.h

  Code template for C functions 
  nearest neighbour of each grid point

  THE FOLLOWING CODE ASSUMES THAT POINT PATTERN (xp, yp) IS SORTED
  IN ASCENDING ORDER OF x COORDINATE

  This code is #included multiple times in nngrid.c 
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
  Either or both DIST and WHICH may be defined.

  Copyright (C) Adrian Baddeley, Jens Oehlschlagel and Rolf Turner 2000-2013
  Licence: GPL >= 2

  $Revision: 1.4 $  $Date: 2014/02/18 08:43:29 $


*/
#endif

void FNAME(nx, x0, xstep,  
	   ny, y0, ystep,   /* pixel grid dimensions */
           np, xp, yp,   /* data points */
	   nnd, nnwhich, 
	   huge)
     /* inputs */
     int *nx, *ny, *np;
     double *x0, *xstep, *y0, *ystep, *huge;
     double *xp, *yp;
     /* outputs */
     double *nnd;
     int *nnwhich;
     /* some inputs + outputs are not used in all functions */
{ 
  int Nxcol, Nyrow, Npoints;
  int i, j, ijpos;
  int mleft, mright, mwhich, lastmwhich;
  double  X0, Y0, Xstep, Ystep;
  double d2, d2min, xj, yi, dx, dy, dx2, hu, hu2;

  Nxcol   = *nx;
  Nyrow   = *ny;
  Npoints = *np;
  hu      = *huge;
  X0      = *x0;
  Y0      = *y0;
  Xstep   = *xstep;
  Ystep   = *ystep;

  hu2      = hu * hu;

  if(Npoints == 0)
    return;

  lastmwhich = 0;

  /* loop over pixels */

  for(j = 0, xj = X0; j < Nxcol; j++, xj += Xstep) {

    R_CheckUserInterrupt();
    
    for(i = 0, yi = Y0; i < Nyrow; i++, yi += Ystep) {

      /* reset nn distance and index */
      d2min = hu2;
      mwhich = -1;

      if(lastmwhich < Npoints) {
	/* search forward from previous nearest neighbour  */
	for(mright = lastmwhich; mright < Npoints; ++mright)
	  {
	    dx = xp[mright] - xj;
	    dx2 = dx * dx; 
	    if(dx2 > d2min) /* note that dx2 >= d2min could break too early */
	      break;
	    dy = yp[mright] - yi;
	    d2 =  dy * dy + dx2;
	    if (d2 < d2min) {
	      /* save as nearest neighbour */
	      d2min = d2;
	      mwhich = mright;
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
	    if(dx2 > d2min) /* note that dx2 >= d2min could break too early */
	      break;
	    dy = yp[mleft] - yi;
	    d2 =  dy * dy + dx2;
	    if (d2 < d2min) {
	      /* save as nearest neighbour */
	      d2min = d2;
	      mwhich = mleft;
	    }
	  }
	/* end backward search */
      }
      /* remember index of most recently-encountered neighbour */
      lastmwhich = mwhich;
      /* copy nn distance for grid point (i, j)
	 to output array nnd[i, j] 
      */
      ijpos = i + j * Nyrow;
#ifdef DIST
      nnd[ijpos] = sqrt(d2min);
#endif
#ifdef WHICH
      nnwhich[ijpos] = mwhich + 1;  /* R indexing */
#endif
    /* end of loop over grid points (i, j) */
    }
  }
}



