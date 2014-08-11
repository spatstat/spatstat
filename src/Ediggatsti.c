#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>
#include "chunkloop.h"

/*

  Ediggatsti.c

  $Revision: 1.2 $     $Date: 2012/03/28 05:55:38 $

  C implementation of 'eval' for DiggleGatesStibbard interaction 

  Assumes point patterns are sorted in increasing order of x coordinate

*/

void Ediggatsti(nnsource, xsource, ysource, idsource, 
		nntarget, xtarget, ytarget, idtarget, 
		rrho, values) 
     /* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget;
     int *idsource, *idtarget;
     double *rrho;
     /* output */
     double *values;
{
  int nsource, ntarget, maxchunk, j, i, ileft, idsourcej;
  double xsourcej, ysourcej, xleft, dx, dy, dx2, d2;
  double rho, rho2, coef, product;

  nsource = *nnsource;
  ntarget = *nntarget;
  rho   = *rrho;

  rho2   = rho * rho;
  coef   = M_PI_2/rho;

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
	 adjust starting position

      */
      xleft  = xsourcej - rho;
      while((xtarget[ileft] < xleft) && (ileft+1 < ntarget))
	++ileft;
      /* 
	 process from ileft until dx > rho
      */
      for(i=ileft; i < ntarget; i++) {
	dx = xtarget[i] - xsourcej;
	dx2 = dx * dx;
	if(dx2 > rho2)
	  break;
	if(idtarget[i] != idsourcej) {
	  dy = ytarget[i] - ysourcej;
	  d2 = dx2 + dy * dy;
	  if(d2 <= rho2) 
	    product *= sin(sqrt(d2) * coef);
	}
      }
      values[j] = log(product * product);
    }
  }
}


