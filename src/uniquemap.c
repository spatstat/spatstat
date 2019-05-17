/*

  uniquemap.c

  !! Assumes points are ordered by increasing x value !!

  $Revision: 1.1 $ $Date: 2019/05/17 07:46:48 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
  Licence: GNU Public Licence >= 2

*/


#include <R.h>

#undef ZCOORD

#undef QUITANY

#define FUNNAME uniqmapxy
#include "uniquemap.h"
#undef FUNNAME

#define QUITANY

#define FUNNAME anydupxy
#include "uniquemap.h"
#undef FUNNAME


