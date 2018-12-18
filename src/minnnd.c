/*

  minnnd.c

  Minimum/Maximum Nearest Neighbour Distance

  Uses code templates in minnnd.h, maxnnd.h

  $Revision: 1.5 $     $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#undef IGNOREZERO

#define FNAME minnnd2
#include "minnnd.h"
#undef FNAME

#define FNAME maxnnd2
#include "maxnnd.h"
#undef FNAME

/* min/max nearest neighbour distance ignoring zero distances */

#define IGNOREZERO

#define FNAME minPnnd2
#include "minnnd.h"
#undef FNAME

#define FNAME maxPnnd2
#include "maxnnd.h"
#undef FNAME


