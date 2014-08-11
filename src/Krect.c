/*

  Krect.c

  $Revision: 1.3 $     $Date: 2014/02/09 03:02:42 $

  +++  Copyright (C) Adrian Baddeley, Julian Gilbey and Rolf Turner 2014 ++++

  Fast code for K function in rectangular case.

     **Assumes point pattern is sorted in increasing order of x coordinate**
     **Assumes window is (0,wide) x (0, high) **
     **Assumes output vectors were initialised to zero**

  Krect.c          defines three interface functions,
                   for weighted, unweighted double, and unweighted integer cases

  KrectFunDec.h    (#included thrice)
                   Function declaration, arguments, storage allocation
  
  KrectV1.h        split according to whether Isotropic Correction is wanted
                   Macro ISOTROPIC is #defined 

  KrectV2.h        split according to whether Translation Correction is wanted
                   Macro TRANSLATION is #defined 

  KrectV3.h        split according to whether Border Correction is wanted
                   Macro BORDER is #defined 

  KrectV4.h        split according to whether Uncorrected estimate is wanted
                   Macro UNCORRECTED is #defined 

  KrectBody.h      Function body, including loops over i and j

  KrectIncrem.h    (#included twice)
                   Code performed when a close pair of points has
                   been found: calculate edge corrections, increment results.

*/

#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

/* This constant is defined in Rmath.h */
#define TWOPI M_2PI

#define ABS(X) (((X) >= 0) ? (X) : (-X))
#define SMALL(X) ((ABS(X) < 1.0e-12) ? 1 : 0)
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))

#undef FNAME
#undef WEIGHTED
#undef COUNTTYPE

#define FNAME KrectInt
#define COUNTTYPE int
#include "KrectFunDec.h"

#undef FNAME
#undef WEIGHTED
#undef COUNTTYPE

#define FNAME KrectDbl
#define COUNTTYPE double
#include "KrectFunDec.h"

#undef FNAME
#undef WEIGHTED
#undef COUNTTYPE

#define FNAME KrectWtd
#define COUNTTYPE double
#define WEIGHTED
#include "KrectFunDec.h"



