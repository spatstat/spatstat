/*
  Function declarations from mhsnoop.c

  $Revision: 1.5 $ $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#include "mhsnoopdef.h"

void initmhsnoop(Snoop *s, SEXP env);

void mhsnoop(Snoop *s, 
	     int irep, 
	     Algor *algo,
	     State *state,
	     Propo *prop,
	     double numer, 
	     double denom,
	     int *itype);

