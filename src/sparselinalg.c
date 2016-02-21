#include <R.h>
#include <R_ext/Utils.h>

/*
  sparselinalg.c

  Counterpart of 'linalg.c' for sparse matrices/arrays

  $Revision: 1.6 $  $Date: 2016/02/20 11:14:12 $

 */

#undef DBG

#define FNAME CspaSumSymOut
#undef WEIGHTS
#include "spasumsymout.h"
#undef FNAME

#define FNAME CspaWtSumSymOut
#define WEIGHTS
#include "spasumsymout.h"
#undef FNAME

