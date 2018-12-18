/*

  minnnd.h

  Code template for minnnd
   to be #included in minnnd.c

  Macros:  
  FNAME          Function name 
  IGNOREZERO     #defined if zero distances should be ignored

  $Revision: 1.2 $     $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/


/* THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER */

void FNAME(n, x, y, huge, result) 
     /* inputs */
     int *n;
     double *x, *y, *huge;
     /* outputs */
     double *result;
{ 
  int npoints, i, maxchunk, left, right;
  double d2, d2min, xi, yi, dx, dy, dy2, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  d2min = hu2;

  if(npoints == 0) return;

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 

  while(i < npoints) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > npoints) maxchunk = npoints;

    for(; i < maxchunk; i++) {

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
#ifdef IGNOREZERO
	      if(d2 > 0) d2min = d2;
#else
	      d2min = d2;
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
#ifdef IGNOREZERO
	      if(d2 > 0) d2min = d2;
#else
	      d2min = d2;
#endif
	  }
	}
      }
    }
  }
  *result = d2min;
}

