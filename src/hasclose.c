/*

  hasclose.c

  $Revision: 1.4 $ $Date: 2016/11/29 05:09:25 $

  Determine whether a point has a neighbour closer than 'r'
  
  Data must be ordered by increasing x coordinate
*/

#include <R.h>

#undef BUG

#undef TORUS

#undef ZCOORD

#define CLOSEFUN hasXclose
#define CROSSFUN hasXYclose
#include "hasclose.h"
#undef CLOSEFUN
#undef CROSSFUN

#define ZCOORD

#define CLOSEFUN hasX3close
#define CROSSFUN hasXY3close
#include "hasclose.h"
#undef CLOSEFUN
#undef CROSSFUN

#define TORUS

#undef ZCOORD

#define CLOSEFUN hasXpclose
#define CROSSFUN hasXYpclose
#include "hasclose.h"
#undef CLOSEFUN
#undef CROSSFUN

#define ZCOORD

#define CLOSEFUN hasX3pclose
#define CROSSFUN hasXY3pclose
#include "hasclose.h"
#undef CLOSEFUN
#undef CROSSFUN
