/*

  nngrid.c

  Nearest Neighbour Distances from a pixel grid to a point pattern

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2013
  Licence: GNU Public Licence >= 2

  $Revision: 1.1 $     $Date: 2013/09/29 08:52:42 $

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

