/* 
   KrectV4.h

   with or without border correction

 */

if((*doBord) == 1) {

#define BORDER
#include "KrectV4.h"

 } else {

#undef BORDER
#include "KrectV4.h"

 }

