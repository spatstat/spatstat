/*

  uniquemap.c

  !! Assumes points are ordered by increasing x value !!

  $Revision: 1.2 $ $Date: 2019/05/21 07:36:34 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
  Licence: GNU Public Licence >= 2

*/


#include <R.h>

#undef ZCOORD

#undef MARKED

#undef QUITANY
#define FUNNAME uniqmapxy
#include "uniquemap.h"
#undef FUNNAME

#define QUITANY
#define FUNNAME anydupxy
#include "uniquemap.h"
#undef FUNNAME
#undef QUITANY

#define MARKED

#undef QUITANY
#define FUNNAME uniqmap2M
#include "uniquemap.h"
#undef FUNNAME

#define QUITANY
#define FUNNAME anydup2M
#include "uniquemap.h"
#undef FUNNAME
#undef QUITANY
