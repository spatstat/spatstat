/* 
   mhv3.h

   tracking or not

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#undef MH_TRACKING

if(tracking) {
  /* saving transition history */
#define MH_TRACKING YES
#include "mhv4.h"
#undef MH_TRACKING
} else {
  /* not saving transition history */
#define MH_TRACKING NO
#include "mhv4.h"
#undef MH_TRACKING
} 
