/*
  inxyp.c

  Point-in-polygon test

  NB: relative to other versions, 'score' is multiplied by 2
  (and is an integer)

  $Revision: 1.7 $   $Date: 2013/09/18 04:20:13 $

 */

#include <R_ext/Utils.h>
#include "chunkloop.h"

void inxyp(x,y,xp,yp,npts,nedges,score,onbndry) 
  /* inputs */
  double *x, *y; /* points to be tested */
  int *npts;    
  double *xp, *yp; /* polygon vertices */
  int *nedges;
  /* outputs */
  int *score;
  int *onbndry;
{
  int i, j, Npts, Nedges, Ne1, contrib, maxchunk;
  double x0, y0, x1, y1, dx, dy, xj, yj, xcrit, ycrit;
  
  Npts = *npts;
  Nedges = *nedges;
  Ne1 = Nedges - 1;

  x0 = xp[Ne1];
  y0 = yp[Ne1];

  OUTERCHUNKLOOP(i, Nedges, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Nedges, maxchunk, 16384) {
      /* visit edge (x0,y0) -> (x1,y1) */
      x1 = xp[i];
      y1 = yp[i];
      dx = x1 - x0;
      dy = y1 - y0;
      for(j = 0; j < Npts; j++) {
	xj = x[j];
	yj = y[j];
	xcrit = (xj - x0) * (xj - x1);
	if(xcrit <= 0) {
	  if(xcrit == 0) {
	    contrib = 1;
	  } else {
	    contrib = 2;
	  }
	  ycrit = yj * dx - xj * dy + x0 * dy - y0 * dx;
	  if(dx < 0) {
	    if(ycrit >= 0)
	      score[j] +=  contrib;
	    onbndry[j] = onbndry[j] | (ycrit == 0);
	  } else if(dx > 0) {
	    if(ycrit < 0) 
	      score[j] -= contrib;
	    onbndry[j] = onbndry[j] | (ycrit == 0);
	  } else {
	    if(xj == x0) 
	      ycrit = (yj - y0) * (yj - y1);
	    onbndry[j] = onbndry[j] | (ycrit <= 0);
	  }
	}
      }
      /* next edge */
      x0 = x1;
      y0 = y1;
    }
  }
}
