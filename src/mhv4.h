/* 
   mhv4.h

   visual debugger or not

*/

#undef MH_SNOOP

if(snooper.active) {
  /* visual debugger */
#define MH_SNOOP YES
#include "mhv5.h"
#undef MH_SNOOP
} else {
  /* no visual debugger */
#define MH_SNOOP NO
#include "mhv5.h"
#undef MH_SNOOP
} 

