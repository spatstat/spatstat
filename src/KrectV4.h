/* 
   KrectV5.h

   with or without uncorrected estimator

 */

if((*doUnco) == 1) {

#define UNCORRECTED
#include "KrectBody.h"

 } else {

#undef UNCORRECTED
#include "KrectBody.h"

 }

