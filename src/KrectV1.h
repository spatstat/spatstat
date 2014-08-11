/* 
   KrectV2.h

   with or without isotropic correction

 */

if((*doIso) == 1) {

#define ISOTROPIC
#include "KrectV2.h"

 } else {

#undef ISOTROPIC
#include "KrectV2.h"

 }

