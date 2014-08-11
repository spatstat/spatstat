/*

  knngrid.c

  K-th Nearest Neighbour Distances from a pixel grid to a point pattern

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2013
  Licence: GNU Public Licence >= 2

  $Revision: 1.6 $     $Date: 2013/11/03 05:06:28 $

  Function body definition is #included from knngrid.h 

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

void knnGinterface(nx, x0, xstep,  
		   ny, y0, ystep,   /* pixel grid dimensions */
		   np, xp, yp,   /* data points */
		   kmax,
		   wantdist, wantwhich,
		   nnd, nnwhich, 
		   huge)
     /* inputs */
     int *nx, *ny, *np;
     double *x0, *xstep, *y0, *ystep, *huge;
     double *xp, *yp;
     int *kmax;
     /* options */
     int *wantdist, *wantwhich;
     /* outputs */
     double *nnd;
     int *nnwhich;
     /* some inputs + outputs are not used in all functions */
{
  void knnGdw(), knnGd(), knnGw();
  int di, wh;
  di = (*wantdist != 0);
  wh = (*wantwhich != 0);
  if(di && wh) {
    knnGdw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, kmax, nnd, nnwhich, huge);
  } else if(di) {
    knnGd(nx, x0, xstep, ny, y0, ystep, np, xp, yp, kmax, nnd, nnwhich, huge);
  } else if(wh) {
    knnGw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, kmax, nnd, nnwhich, huge);
  }
}

#undef FNAME
#undef DIST
#undef WHICH

/* 
   knnGdw

   nearest neighbours 1:kmax

   returns distances and indices

*/

#define FNAME knnGdw
#define DIST
#define WHICH
#include "knngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* 
   knnGd

   nearest neighbours 1:kmax

   returns distances only

*/

#define FNAME knnGd
#define DIST
#include "knngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* 
   knnGw 

   nearest neighbours 1:kmax

   returns indices only

*/

#define FNAME knnGw
#define WHICH
#include "knngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

