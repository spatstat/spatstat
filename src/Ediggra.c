#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/*

  Ediggra.c

  $Revision: 1.4 $     $Date: 2012/03/28 05:55:50 $

  C implementation of 'eval' for DiggleGratton interaction (exponentiated)

  Assumes point patterns are sorted in increasing order of x coordinate

*/

double sqrt();

void Ediggra(nnsource, xsource, ysource, idsource, 
	     nntarget, xtarget, ytarget, idtarget, 
	     ddelta, rrho, values) 
     /* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget;
     int *idsource, *idtarget;
     double *ddelta, *rrho;
     /* output */
     double *values;
{
  int nsource, ntarget, maxchunk, j, i, ileft, idsourcej;
  double xsourcej, ysourcej, xleft, dx, dy, dx2, d2;
  double delta, rho, delta2, rho2, rhominusdelta;
  double product;

  nsource = *nnsource;
  ntarget = *nntarget;
  delta = *ddelta;
  rho   = *rrho;

  rho2   = rho * rho;
  delta2 = delta * delta;
  rhominusdelta = rho - delta;

  if(nsource == 0 || ntarget == 0) 
    return;

  ileft = 0;

  OUTERCHUNKLOOP(j, nsource, maxchunk, 65536) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nsource, maxchunk, 65536) {
      product = 1;
      xsourcej = xsource[j];
      ysourcej = ysource[j];
      idsourcej = idsource[j];

      /* 
	 adjust starting point
      */

      xleft  = xsourcej - rho;
      while((xtarget[ileft] < xleft) && (ileft+1 < ntarget))
	++ileft;

      /* 
	 process until dx > rho 
	 (or until product is zero)
      */
      for(i=ileft; i < ntarget; i++) {
	dx = xtarget[i] - xsourcej;
	dx2 = dx * dx;
	if(dx2 > rho2) 
	  break;
	if(idtarget[i] != idsourcej) {
	  dy = ytarget[i] - ysourcej;
	  d2 = dx2 + dy * dy;
	  if(d2 <= rho2) {
	    if(d2 <= delta2) {
	      product = 0;
	      break;
	    }
	    else 
	      product *= (sqrt(d2) - delta)/rhominusdelta;
	  }
	}
      }
      values[j] = product;
    }
  }
}
