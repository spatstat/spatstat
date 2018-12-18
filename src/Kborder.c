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
   
  $Revision: 1.5 $ $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

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



