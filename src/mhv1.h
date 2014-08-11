/* 
   mhv1.h

   marked or unmarked simulation

*/

#undef MH_MARKED

if(marked) {
  /* marked process */
#define MH_MARKED TRUE
#include "mhv2.h"
#undef MH_MARKED
} else {
  /* unmarked process */
#define MH_MARKED FALSE
#include "mhv2.h"
#undef MH_MARKED
} 
