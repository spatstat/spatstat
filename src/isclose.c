/*

  isclose.c

  $Revision: 1.1 $ $Date: 2016/10/14 01:22:03 $

  Determine whether a point has a neighbour closer than 'r'
  
  Data must be ordered by increasing x coordinate
*/

#include <R.h>

#define CLOSEFUN isXclose
#define CROSSFUN isXYclose
#undef ZCOORD
#include "isclose.h"
#undef CLOSEFUN
#undef CROSSFUN

#define CLOSEFUN isX3close
#define CROSSFUN isXY3close
#define ZCOORD
#include "isclose.h"
#undef CLOSEFUN
#undef CROSSFUN
#undef ZCOORD
