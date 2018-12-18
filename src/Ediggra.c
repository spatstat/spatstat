#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"
#include "looptest.h"

/*

  Ediggra.c

  $Revision: 1.8 $     $Date: 2018/12/18 02:43:11 $

  C implementation of 'eval' for DiggleGratton interaction (exponentiated)

  Assumes point patterns are sorted in increasing order of x coordinate

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

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
  double delta, rho, delta2, rho2, rho2pluseps, rhominusdelta;
  double product;

  nsource = *nnsource;
  ntarget = *nntarget;
  delta = *ddelta;
  rho   = *rrho;

  if(nsource == 0 || ntarget == 0) 
    return;

  rho2   = rho * rho;
  delta2 = delta * delta;
  rhominusdelta = rho - delta;
  rho2pluseps = rho2 + EPSILON(rho2);

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
	if(dx2 > rho2pluseps) 
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


/* 'split' version separating hard core terms from others */

void ESdiggra(nnsource, xsource, ysource, idsource, 
	      nntarget, xtarget, ytarget, idtarget, 
	      ddelta, rrho, positive, hardcore) 
     /* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget;
     int *idsource, *idtarget;
     double *ddelta, *rrho;
     /* output */
     double *positive;
     int *hardcore;
{
  int nsource, ntarget, maxchunk, j, i, ileft, idsourcej;
  double xsourcej, ysourcej, xleft, dx, dy, dx2, d2;
  double delta, rho, delta2, rho2, rho2pluseps, rhominusdelta;
  double product;

  nsource = *nnsource;
  ntarget = *nntarget;
  delta = *ddelta;
  rho   = *rrho;

  if(nsource == 0 || ntarget == 0) 
    return;

  rho2   = rho * rho;
  delta2 = delta * delta;
  rhominusdelta = rho - delta;
  rho2pluseps = rho2 + EPSILON(rho2);

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
      */
      for(i=ileft; i < ntarget; i++) {
	dx = xtarget[i] - xsourcej;
	dx2 = dx * dx;
	if(dx2 > rho2pluseps) 
	  break;
	if(idtarget[i] != idsourcej) {
	  dy = ytarget[i] - ysourcej;
	  d2 = dx2 + dy * dy;
	  if(d2 <= rho2) {
	    if(d2 <= delta2) {
	      hardcore[j] = 1;
	    } else {
	      product *= (sqrt(d2) - delta)/rhominusdelta;
	    }
	  }
	}
      }
      positive[j] = product;
    }
  }
}

