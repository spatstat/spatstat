/*
  nndist.h

  Code template for C functions supporting nndist and nnwhich (k=1)

  THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER

  This code is #included multiple times in nndistance.c 
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
  Either or both DIST and WHICH may be defined.

  Copyright (C) Adrian Baddeley, Jens Oehlschlagel and Rolf Turner 2000-2012
  Licence: GPL >= 2

  $Revision: 1.2 $  $Date: 2012/03/14 02:37:27 $

*/

void FNAME(n, x, y, 
#ifdef DIST
	   nnd,
#endif
#ifdef WHICH
	   nnwhich, 
#endif
           huge)
     /* inputs */
     int *n;
     double *x, *y, *huge;
     /* outputs */
#ifdef DIST
     double *nnd;
#endif
#ifdef WHICH
     int *nnwhich;
#endif
{ 
  int npoints, i, maxchunk, left, right;
  double d2, d2min, xi, yi, dx, dy, dy2, hu, hu2;
#ifdef WHICH
  int which;
#endif

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < npoints) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > npoints) maxchunk = npoints;

    for(; i < maxchunk; i++) {

      d2min = hu2;
#ifdef WHICH
      which = -1;
#endif
      xi = x[i];
      yi = y[i];

      if(i < npoints - 1) {
	/* search forward */
	for(right = i + 1; right < npoints; ++right)
	  {
	    dy = y[right] - yi;
	    dy2 = dy * dy;
	    if(dy2 > d2min)
	      break;
	    dx = x[right] - xi;
	    d2 =  dx * dx + dy2;
	    if (d2 < d2min) {
	      d2min = d2;
#ifdef WHICH
	      which = right;
#endif
	    }
	  }
      }
      if(i > 0){
	/* search backward */
	for(left = i - 1; left >= 0; --left)
	{
	  dy = yi - y[left];
	  dy2 = dy * dy;
	  if(dy2 > d2min)
	    break;

	  dx = x[left] - xi;
	  d2 =  dx * dx + dy2;
	  if (d2 < d2min) {
	    d2min = d2;
#ifdef WHICH
	    which = left;
#endif
	  }
	}
      }

#ifdef DIST
      nnd[i] = sqrt(d2min);
#endif
#ifdef WHICH
      nnwhich[i] = which + 1; /* R indexing */
#endif
    }
  }
}
