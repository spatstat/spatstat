#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"

/*

  localpcf.c

  $Revision: 1.4 $     $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  Assumes point patterns are sorted in increasing order of x coordinate

*/

#undef WEIGHTED

#include "localpcf.h"

#define WEIGHTED 1

#include "localpcf.h"
