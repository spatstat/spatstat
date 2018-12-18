#include <R.h>
#include "yesno.h"

/*

  linknnd.c

  k-th nearest neighbours in a linear network

  Sparse representation of network

  ! Data points must be ordered by segment index !


  $Revision: 1.4 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

#undef HUH

#undef CROSS
#define FNAME linknnd
#include "linknnd.h"
#undef FNAME

#define CROSS
#define FNAME linknncross
#include "linknnd.h"
#undef CROSS
#undef FNAME

