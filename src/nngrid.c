/*

  nngrid.c

  Nearest Neighbour Distances from a pixel grid to a point pattern

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2013
  Licence: GNU Public Licence >= 2

  $Revision: 1.4 $     $Date: 2013/11/03 03:41:23 $

  Function body definition is #included from nngrid.h 

  THE FOLLOWING FUNCTIONS ASSUME THAT x IS SORTED IN ASCENDING ORDER 

*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "yesno.h"

double sqrt();

/* THE FOLLOWING CODE ASSUMES THAT x IS SORTED IN ASCENDING ORDER */

/* general interface */

void nnGinterface(nx, x0, xstep,  
		  ny, y0, ystep,   /* pixel grid dimensions */
		  np, xp, yp,   /* data points */
		  wantdist, wantwhich, /* options */
		  nnd, nnwhich, 
		  huge)
     /* inputs */
     int *nx, *ny, *np;
     double *x0, *xstep, *y0, *ystep, *huge;
     double *xp, *yp;
     /* options */
     int *wantdist, *wantwhich;
     /* outputs */
     double *nnd;
     int *nnwhich;
     /* some inputs + outputs are not used in all functions */
{
  void nnGdw(), nnGd(), nnGw();
  int di, wh;
  di = (*wantdist != 0);
  wh = (*wantwhich != 0);
  if(di && wh) {
    nnGdw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
  } else if(di) {
    nnGd(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
  } else if(wh) {
    nnGw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
  }
}


#undef FNAME
#undef DIST
#undef WHICH

/* 
   nnGdw

   returns distances and indices

*/

#define FNAME nnGdw
#define DIST
#define WHICH
#include "nngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* 
   nnGd

   returns distances only

*/

#define FNAME nnGd
#define DIST
#include "nngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* 
   nnGw 

   returns indices only

*/

#define FNAME nnGw
#define WHICH
#include "nngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

