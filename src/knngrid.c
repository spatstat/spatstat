/*

  knngrid.c

  K-th Nearest Neighbour Distances from a pixel grid to a point pattern

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2013
  Licence: GNU Public Licence >= 2

  $Revision: 1.2 $     $Date: 2013/05/27 02:09:10 $

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

