/*

  close3pair.c

  $Revision: 1.31 $     $Date: 2013/11/22 01:01:12 $

  closepairs and crosspairs for 3D

  Assumes point pattern is sorted in increasing order of x coordinate

*/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

#define OK 0
#define ERR_OVERFLOW 1
#define ERR_ALLOC 2

#define intRealloc(PTR, OLDLENGTH, NEWLENGTH) \
  (int *) S_realloc((char *) PTR, NEWLENGTH, OLDLENGTH, sizeof(int))

#define dblRealloc(PTR, OLDLENGTH, NEWLENGTH) \
  (double *) S_realloc((char *) PTR, NEWLENGTH, OLDLENGTH, sizeof(double))

double sqrt();

/* .......  define functions, using closefuns.h  ........*/

/* enable 3D code */
#define ZCOORD

/* return i, j only */
#define CLOSEFUN close3IJpairs
#define CROSSFUN cross3IJpairs
#undef THRESH
#undef COORDS
#include "closefuns.h"
#undef CLOSEFUN
#undef CROSSFUN
#undef THRESH
#undef COORDS

/* return i, j, xi, yi, zi, xj, yj, zj, dx, dy, dz, d */
#define CLOSEFUN close3pairs
#define CROSSFUN cross3pairs
#undef THRESH
#define COORDS
#include "closefuns.h"
#undef CLOSEFUN
#undef CROSSFUN
#undef THRESH
#undef COORDS

/* return i, j, t where t = 1{d < s} */

#define CLOSEFUN close3thresh
#define CROSSFUN cross3thresh
#define THRESH
#undef COORDS
#include "closefuns.h"
#undef CLOSEFUN
#undef CROSSFUN
#undef THRESH
#undef COORDS

