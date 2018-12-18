#include <R.h>
#include "yesno.h"

/*

  linvdist.c

  Distance function at vertices
  (shortest distance from each vertex to a data point)

  Sparse representation of network

  $Revision: 1.1 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  ! Data points must be ordered by segment index !

 */

#undef HUH

/* definition of Clinvdist */
#define FNAME Clinvdist
#undef WHICH
#include "linvdist.h"

/* definition of Clinvwhichdist */
#undef FNAME
#define FNAME Clinvwhichdist
#define WHICH
#include "linvdist.h"
