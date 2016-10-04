#include <R.h>
#include <math.h>
#include "yesno.h"

/* 
   linequad.c

   make a quadrature scheme on a linear network

   Clinequad    unmarked pattern
   ClineMquad   multitype pattern

   $Revision: 1.5 $  $Date: 2016/10/03 08:43:57 $ 

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

