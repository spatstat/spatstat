/*

  knndistance.c

  K-th Nearest Neighbour Distances between points

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2013
  Licence: GNU Public Licence >= 2

  $Revision: 1.3 $     $Date: 2013/05/27 02:09:10 $

  Function definitions are #included from knndist.h and knnXdist.h

  THE FOLLOWING FUNCTIONS ASSUME THAT y IS SORTED IN ASCENDING ORDER 

  SINGLE LIST:
  knndsort     k-th nearest neighbour distances
  knnwhich     k-th nearest neighbours
  knnsort      k-th nearest neighbours and their distances

  ONE LIST TO ANOTHER LIST:
  knnXdist     Nearest neighbour distance from one list to another
  knnXwhich    Nearest neighbour ID from one list to another
  knnX         Nearest neighbour ID & distance from one list to another

  ONE LIST TO ANOTHER OVERLAPPING LIST:
  knnXEdist    Nearest neighbour distance from one list to another, overlapping
  knnXEwhich   Nearest neighbour ID from one list to another, overlapping
  knnXE        Nearest neighbour ID & distance 


*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "yesno.h"

double sqrt();

/* THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER */

/* ------------------- one point pattern X --------------------- */

/* 
   knndsort 

   nearest neighbours 1:kmax

   returns distances only

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knndsort
#define DIST
#include "knndist.h"

/* 
   knnwhich

   nearest neighbours 1:kmax

   returns identifiers only

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knnwhich
#define WHICH
#include "knndist.h"

/* 
   knnsort 

   nearest neighbours 1:kmax

   returns distances and indices

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knnsort
#define DIST
#define WHICH
#include "knndist.h"

/* --------------- two distinct point patterns X and Y --------------- */

/* 
   knnXdist

   returns distances only

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knnXdist
#define DIST
#include "knnXdist.h"

/* 
   knnXwhich

   returns identifiers only

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knnXwhich
#define WHICH
#include "knnXdist.h"

/* 
   knnX 

   returns distances and indices

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knnX
#define DIST
#define WHICH
#include "knnXdist.h"

/* --------------- overlapping point patterns X and Y --------------- */

/* 
   knnXEdist

   returns distances only

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knnXEdist
#define DIST
#define EXCLUDE
#include "knnXdist.h"

/* 
   knnXEwhich

   returns identifiers only

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knnXEwhich
#define WHICH
#define EXCLUDE
#include "knnXdist.h"

/* 
   knnXE 

   returns distances and indices

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knnXE
#define DIST
#define WHICH
#define EXCLUDE
#include "knnXdist.h"

