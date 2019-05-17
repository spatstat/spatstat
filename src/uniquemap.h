/*
  uniquemap.h

  Function definitions to be #included in uniquemap.c
  several times with different values of macros.

  !! Assumes points are ordered by increasing x value !!

  Assumes <R.h> is included

  Macros used:

  FUNNAME    name of function

  QUITANY    return TRUE immediately if any duplicates are found

  ZCOORD     if defined, coordinates are 3-dimensional

  $Revision: 1.5 $ $Date: 2019/05/17 08:19:25 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
  Licence: GNU Public Licence >= 2

*/

#ifdef ZCOORD 
#define SPACEDIM 3
#else
#define SPACEDIM 2
#endif

void FUNNAME(int *n,
	     double *x,
	     double *y,
#ifdef ZCOORD
	     double *z,
#endif
#ifdef QUITANY
	     int *anydup
#else
	     int *uniqmap
#endif
	     )
{
  double xi, yi, dx, dy, d2;
#ifdef ZCOORD
  double zi, dz;
#endif

  int N, maxchunk, i, j;

  /* loop in chunks of 2^16 */

  N = *n;
  
  i = 0; maxchunk = 0; 
  while(i < N) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > N) maxchunk = N;

    for(; i < maxchunk; i++) {

      xi = x[i];
      yi = y[i];
#ifdef ZCOORD
      zi = z[i];
#endif

      if(i + 1 < N) {
#ifndef QUITANY
	if(uniqmap[i] == 0) { /* i.e. don't seek duplicates of a duplicate */
#endif	  
	  /* scan forward */
	  for(j = i + 1; j < N; j++) {
	    dx = x[j] - xi;
	    if(dx > DOUBLE_EPS) 
	      break;
	    dy = y[j] - yi;
	    d2 = dx * dx + dy * dy;
#ifdef ZCOORD
	    if(d2 <= 0.0) {
	      dz = z[j] - zi;
	      d2 = d2 + dz * dz;
#endif
	      if(d2 <= 0.0) {
		/* j is a duplicate of i */
#ifdef QUITANY
		*anydup = 1; /* i.e. TRUE */
		return;
#else	      
		uniqmap[j] = i + 1; /* R indexing */
#endif	      
	      }
#ifdef ZCOORD
	    }
#endif
	  }
#ifndef QUITANY	  
	}
#endif	
      }
    }
  }
}


