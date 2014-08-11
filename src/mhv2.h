/* 
   mhv2.h

   single interaction or hybrid

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

