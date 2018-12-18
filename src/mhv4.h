/* 
   mhv4.h

   visual debugger or not

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

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

