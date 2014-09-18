/*

  minnnd.c

  Minimum/Maximum Nearest Neighbour Distance

  Uses code templates in minnnd.h, maxnnd.h

  $Revision: 1.4 $     $Date: 2014/09/18 01:28:48 $

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


