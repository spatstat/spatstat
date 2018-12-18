/*

   mhsnoopdef.h

   Define structure 'Snoop' containing visual debugger parameters and state

   $Revision: 1.3 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#ifndef R_INTERNALS_H_
#include <Rinternals.h>
#endif

typedef struct Snoop {
  int active;    /* true or false */
  int nextstop;  /* jump to iteration number 'nextstop' */
  int nexttype;  /* jump to the next proposal of type 'nexttype' */
  SEXP env;    /* environment for exchanging data with R */
  SEXP expr;   /* callback expression for visual debugger */
} Snoop;

#define NO_TYPE -1
