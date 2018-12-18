/*

  nn3Ddist.h

  Code template for nearest-neighbour algorithms for 3D point patterns

  Input is a single point pattern - supports 'nndist' and 'nnwhich'

  This code is #included multiple times in nn3Ddist.c
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
  Either or both DIST and WHICH may be defined.

  THE FOLLOWING CODE ASSUMES THAT THE POINT PATTERN IS SORTED
  IN ASCENDING ORDER OF THE z COORDINATE

  $Revision: 1.6 $ $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

  
void FNAME(n, x, y, z, 
	   nnd, nnwhich, huge)
/* inputs */
     int *n;
     double *x, *y, *z, *huge;
     /* outputs */
     double *nnd; 
     int *nnwhich;
{ 
  int npoints, i, j, maxchunk;
  double d2, d2min, xi, yi, zi, dx, dy, dz, dz2, hu, hu2;
#ifdef WHICH
  int which;
#endif

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {
      d2min = hu2;
#ifdef WHICH
      which = -1;
#endif
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
#ifdef WHICH
	    which = j;
#endif
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
#ifdef WHICH
	    which = j;
#endif
	  }
	}
      }
#ifdef DIST
      nnd[i] = sqrt(d2min);
#endif
#ifdef WHICH
      /* convert to R indexing */
      nnwhich[i] = which + 1;
#endif
    }
  }
}

