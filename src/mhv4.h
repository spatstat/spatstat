/* 
   mhv4.h

   visual debugger or not

*/

#undef MH_SNOOP

if(snooper.active) {
  /* visual debugger */
#define MH_SNOOP TRUE
#include "mhloop.h"
#undef MH_SNOOP
} else {
  /* no visual debugger */
#define MH_SNOOP FALSE
#include "mhloop.h"
#undef MH_SNOOP
} 

