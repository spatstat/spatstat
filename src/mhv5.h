/* 
   mhv5.h

   tempered or not

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#undef MH_TEMPER

if(tempered) {
  /* tempering */
#define MH_TEMPER YES
#include "mhloop.h"
#undef MH_TEMPER
} else {
  /* usual, no tempering */
#define MH_TEMPER NO
#include "mhloop.h"
#undef MH_TEMPER
} 

