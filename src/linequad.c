#include <R.h>
#include <math.h>
#include "yesno.h"

/* 
   linequad.c

   make a quadrature scheme on a linear network

   Clinequad    unmarked pattern
   ClineMquad   multitype pattern

   $Revision: 1.6 $  $Date: 2018/12/18 02:43:11 $ 

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

#define SWAP(X,Y,TMP) TMP = Y; Y = X; X = TMP

#undef HUH

#define FUNNAME Clinequad
#define FMKNAME ClineMquad
#undef ALEA
#include "linequad.h"
#undef FUNNAME
#undef FMKNAME

#define FUNNAME ClineRquad
#define FMKNAME ClineRMquad
#define ALEA
#include "linequad.h"
#undef FUNNAME
#undef FMKNAME
#undef ALEA

