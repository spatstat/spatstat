#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/*

  Estrauss.c

  $Revision: 1.3 $     $Date: 2012/03/28 05:56:24 $

  C implementation of 'eval' for Strauss interaction

  Calculates number of data points within distance r of each quadrature point
  (when 'source' = quadrature points, 'target' = data points)

  Assumes point patterns are sorted in increasing order of x coordinate

*/

double sqrt();

void closepaircounts(nnsource, xsource, ysource, 
		     nntarget, xtarget, ytarget, 
		     rrmax, counts) 
/* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget, *rrmax;
     /* output */
     int *counts;
{
  int nsource, ntarget, maxchunk, j, i, ileft, counted;
  double xsourcej, ysourcej, rmax, r2max, xleft, dx, dy, dx2, d2;

  nsource = *nnsource;
  ntarget = *nntarget;
  rmax = *rrmax;
  r2max = rmax * rmax;

  if(nsource == 0 || ntarget == 0) 
    return;

  ileft = 0;

  OUTERCHUNKLOOP(j, nsource, maxchunk, 65536) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nsource, maxchunk, 65536) {
      counted = 0;
      xsourcej = xsource[j];
      ysourcej = ysource[j];
      /* 
	 adjust starting point
      */
      xleft  = xsourcej - rmax;
      while((xtarget[ileft] < xleft) && (ileft+1 < ntarget))
	++ileft;

      /* 
	 process from ileft to iright
      */
      for(i=ileft; i < ntarget; i++) {
	dx = xtarget[i] - xsourcej;
	dx2 = dx * dx;
	if(dx2 > r2max)
	  break;
	dy = ytarget[i] - ysourcej;
	d2 = dx2 + dy * dy;
	if(d2 <= r2max)
	  ++counted;
      }
      counts[j] = counted;
    }
  }
}
