/* 
   mhv2.h

   single interaction or hybrid

*/

#undef MH_SINGLE

if(Ncif == 1) {
  /* single interaction */
#define MH_SINGLE TRUE
#include "mhv3.h"
#undef MH_SINGLE
} else {
  /* hybrid interaction */
#define MH_SINGLE FALSE
#include "mhv3.h"
#undef MH_SINGLE
} 

