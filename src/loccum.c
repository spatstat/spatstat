#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/*

  loccum.c

  $Revision: 1.1 $     $Date: 2013/05/27 02:09:10 $

  Compute local cumulative sums or products of weights

  locsum:  f_i(t) = \sum_{j: j \neq i, ||x_j - x_i|| \le t} v(x_j)
            for a data point pattern {x_i} 

  locxsum: f_u(t) = \sum_{||x_i - u|| \le t} v(x_i)
            for a grid of points {u} and a data point pattern {x_i} 
	    (also works if {u} is another point pattern)

  locprod:  f_i(t) = \prod_{j: j \neq i, ||x_j - x_i|| \le t} v(x_j)
            for a data point pattern {x_i} 

  locxprod: f_u(t) = \prod_{||x_i - u|| \le t} v(x_i)
            for a grid of points {u} and a data point pattern {x_i} 
	    (also works if {u} is another point pattern)

  Assumes point patterns are sorted in increasing order of x coordinate

  Uses C code template files : loccums.h, loccumx.h

*/

/* data-to-data */

#undef FNAME
#undef NULVAL
#undef INC

#define FNAME locsum
#define NULVAL 0.0
#define INC(A,B) A += B

#include "loccums.h"

#undef FNAME
#undef NULVAL
#undef INC

#define FNAME locprod
#define NULVAL 1.0
#define INC(A,B) A *= B

#include "loccums.h"

/* test-grid-to-data */

#undef FNAME
#undef NULVAL
#undef INC

#define FNAME locxsum
#define NULVAL 0.0
#define INC(A,B) A += B

#include "loccumx.h"

#undef FNAME
#undef NULVAL
#undef INC

#define FNAME locxprod
#define NULVAL 1.0
#define INC(A,B) A *= B

#include "loccumx.h"



