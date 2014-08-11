/*

  nndistance.c

  Nearest Neighbour Distances between points

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2012
  Licence: GNU Public Licence >= 2

  $Revision: 1.15 $     $Date: 2012/03/18 07:11:23 $

  THE FOLLOWING FUNCTIONS ASSUME THAT y IS SORTED IN ASCENDING ORDER 

  nndistsort    Nearest neighbour distances 
  nnwhichsort   Nearest neighbours
  nnsort        Nearest neighbours & distances

  nnXwhich      Nearest neighbour from one list to another
  nnXexclude    Nearest neighbour from one list to another, with overlaps

  knndsort      k-th nearest neighbour distances
  knnsort       k-th nearest neighbours and their distances
*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

double sqrt();

/* THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER */

/* ------------------- one point pattern X --------------------- */

/* 
   nndistsort: nearest neighbour distances 
*/

#undef FNAME
#undef DIST
#undef WHICH
#define FNAME nndistsort
#define DIST
#include "nndist.h"

/* 
   nnwhichsort: id of nearest neighbour 
*/

#undef FNAME
#undef DIST
#undef WHICH
#define FNAME nnwhichsort
#define WHICH
#include "nndist.h"

/* 
   nnsort: distance & id of nearest neighbour 
*/

#undef FNAME
#undef DIST
#undef WHICH
#define FNAME nnsort
#define DIST
#define WHICH
#include "nndist.h"

/* --------------- two distinct point patterns X and Y  ----------------- */

/* 
   nnXdist:  nearest neighbour distance
	      (from each point of X to the nearest point of Y)
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXdist
#define DIST
#include "nndistX.h"

/* 
   nnXwhich:  nearest neighbour id
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXwhich
#define WHICH
#include "nndistX.h"

/* 
   nnX:  nearest neighbour distance and id
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnX
#define DIST
#define WHICH
#include "nndistX.h"

/* --------------- two point patterns X and Y with common points --------- */

/*
   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].
*/

/* 
   nnXEdist:  similar to nnXdist
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXEdist
#define DIST
#define EXCLUDE
#include "nndistX.h"

/* 
   nnXEwhich:  similar to nnXwhich
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXEwhich
#define WHICH
#define EXCLUDE
#include "nndistX.h"

/* 
   nnXE:  similar to nnX
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXE
#define DIST
#define WHICH
#define EXCLUDE
#include "nndistX.h"


/*    -------------- k-th nearest neighbours ------------------ */

/* 
   knndsort 

   nearest neighbours 1:kmax

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME knndsort
#define DIST
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

