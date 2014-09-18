/*

  maxnnd.h

  Code template for maxnnd 
   to be #included in minnnd.c

  Macros: 
  FNAME          Function name
  IGNOREZERO     #defined if zero distances should be ignored

  $Revision: 1.2 $     $Date: 2014/09/18 01:00:30 $

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
  double d2, d2mini, d2max, xi, yi, dx, dy, dy2, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  /* maximum (over all i) nearest-neighbour distance, squared */
  d2max = 0.0;

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
      
      /* nearest-neighbour distance for point i,   squared */
      d2mini = hu2;

      if(i < npoints - 1) {
	/* search forward */
	for(right = i + 1; right < npoints; ++right)
	  {
	    dy = y[right] - yi;
	    dy2 = dy * dy;
	    if(dy2 > d2mini)
	      break;
	    dx = x[right] - xi;
	    d2 =  dx * dx + dy2;
	    if (d2 < d2mini) {
#ifdef IGNOREZERO
	      if(d2 > 0) {
#endif
		d2mini = d2;
		if(d2mini <= d2max)
		  break;
#ifdef IGNOREZERO
	      }
#endif
	    }
	  }
      }
      if(i > 0 && d2mini > d2max){
	/* search backward */
	for(left = i - 1; left >= 0; --left)
	{
	  dy = yi - y[left];
	  dy2 = dy * dy;
	  if(dy2 > d2mini)
	    break;

	  dx = x[left] - xi;
	  d2 =  dx * dx + dy2;
	  if (d2 < d2mini) {
#ifdef IGNOREZERO
	    if(d2 > 0) {
#endif
	      d2mini = d2;
	      if(d2mini <= d2max)
		break;
#ifdef IGNOREZERO
	    }
#endif
	  }
	}
      }
      if(d2mini > d2max)
	d2max = d2mini;
    }
  }
  *result = d2max;
}
