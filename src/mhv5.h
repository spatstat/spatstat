/* 
   mhv5.h

   tempered or not

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

