/* 
   mhv1.h

   marked or unmarked simulation

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#undef MH_MARKED

if(marked) {
  /* marked process */
#define MH_MARKED YES
#include "mhv2.h"
#undef MH_MARKED
} else {
  /* unmarked process */
#define MH_MARKED NO
#include "mhv2.h"
#undef MH_MARKED
} 
