/* 
   mhv3.h

   tracking or not

*/

#undef MH_TRACKING

if(tracking) {
  /* saving transition history */
#define MH_TRACKING TRUE
#include "mhv4.h"
#undef MH_TRACKING
} else {
  /* not saving transition history */
#define MH_TRACKING FALSE
#include "mhv4.h"
#undef MH_TRACKING
} 
