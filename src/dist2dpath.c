#include <R.h>
#include <R_ext/Utils.h>

/*
  given matrix of edge lengths
  compute matrix of shortest-path distances

  Uses dist2dpath.h
*/

#define FNAME Ddist2dpath
#define DTYPE double
#define FLOATY

#include "dist2dpath.h"

#undef FNAME
#undef DTYPE
#undef FLOATY

#define FNAME Idist2dpath
#define DTYPE int

#include "dist2dpath.h"

