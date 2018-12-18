#include <R.h>
#include <R_ext/Utils.h>

/*
  sparselinalg.c

  Counterpart of 'linalg.c' for sparse matrices/arrays

  $Revision: 1.7 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

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

