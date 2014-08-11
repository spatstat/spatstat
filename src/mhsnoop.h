/*
  Function declarations from mhsnoop.c

  $Revision: 1.4 $ $Date: 2013/05/27 02:09:10 $

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

