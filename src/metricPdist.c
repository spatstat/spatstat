/*
  metricPdist.c

  Distance transform of binary pixel image
  using arbitrary metric

  This code #includes metricPdist.h multiple times.

  $Revision: 1.1 $  $Date: 2018/05/15 08:50:24 $

 */

/* Once-only declarations: */

#include <math.h>
#include "raster.h"

void   dist_to_bdry();
void   shape_raster();

#define UNDEFINED -1
#define Is_Defined(I) (I >= 0)
#define Is_Undefined(I) (I < 0)

/*

  Definitions for each metric

  For each definition we need the following macros:
  FNAME          Name of the function (that will be called from R)
  MARGLIST       List of arguments to FNAME specifying the metric
  MARGDECLARE    Declarations of these function arguments 
  MTEMPDECLARE   Declaration and initialisation of variables for use by metric
  METRIC         Expression for calculating the metric (x1,y1,x2,y2)
*/

/* 
   (1) Rectangular metric

   Unit ball is a rectangle with width 1 unit, height 'aspect' units.

   mdt  = metric distance transform
   P    = pixel image input
   O    = orthogonally oriented to axis
   rect = rectangular 
*/

#define FNAME mdtPOrect
#define MARGLIST aspect
#define MARGDECLARE double *aspect
#define MTEMPDECLARE double asp; asp=*aspect
#define METRIC(X,Y,XX,YY) rectdist(X,Y,XX,YY,asp)

double rectdist(x, y, xx, yy, asp)
     double x, y, xx, yy, asp;
{
  double dx, dy, d;
  dx = x-xx;
  dy = (y-yy)/asp;
  if(dx < 0) dx = -dx;
  if(dy < 0) dy = -dy;
  d = (dx > dy)? dx :  dy;
  return d;
}

#include "metricPdist.h"

#undef FNAME
#undef MARGLIST
#undef MARGDECLARE
#undef MTEMPDECLARE 
#undef METRIC


