/* 
   mhv1.h

   marked or unmarked simulation

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
