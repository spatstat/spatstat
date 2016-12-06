#include <R.h>
#include "yesno.h"

/*

  linknnd.c

  k-th nearest neighbours in a linear network

  Sparse representation of network

  ! Data points must be ordered by segment index !


  $Revision: 1.3 $  $Date: 2016/12/04 11:08:58 $

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

