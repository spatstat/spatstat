/*
  Function declarations from mhsnoop.c

  $Revision: 1.3 $ $Date: 2013/02/12 05:27:46 $

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

