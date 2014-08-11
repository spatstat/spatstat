#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

/* 
  Kborder.c

  Efficient computation of border-corrected estimates of K 
  for large datasets
  
  KborderI()   Estimates K function, 
               returns integer numerator & denominator

  KborderD()   Estimates K function, 
               returns double precision numerator & denominator

  Kwborder()   Estimates Kinhom.

  Functions require (x,y) data to be sorted in ascending order of x
  and expect r values to be equally spaced and starting at zero
   
  $Revision: 1.3 $ $Date: 2012/02/28 10:55:59 $

*/

#undef WEIGHTED

#define FNAME KborderI
#define OUTTYPE int
#include "Kborder.h"

#undef FNAME
#undef OUTTYPE
#define FNAME KborderD
#define OUTTYPE double
#include "Kborder.h"

#undef FNAME
#undef OUTTYPE
#define FNAME Kwborder
#define WEIGHTED
#define OUTTYPE double
#include "Kborder.h"



