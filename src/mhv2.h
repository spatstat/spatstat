/* 
   mhv2.h

   single interaction or hybrid

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#undef MH_SINGLE

if(Ncif == 1) {
  /* single interaction */
#define MH_SINGLE YES
#include "mhv3.h"
#undef MH_SINGLE
} else {
  /* hybrid interaction */
#define MH_SINGLE NO
#include "mhv3.h"
#undef MH_SINGLE
} 

