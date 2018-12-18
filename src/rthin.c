#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

/* 
   rthin.c

   Select from the integers 1:n with probability p
   by simulating geometric(p) jumps between selected integers

   $Revision: 1.2 $ $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

SEXP thinjumpequal(SEXP n,
		   SEXP p,
		   SEXP guess) 
{
  int N;
  double P;

  int *w;  /* temporary storage for selected integers */
  int nw, nwmax;

  int i, j, k;
  double log1u, log1p;

  /* R object return value */
  SEXP Out;
  /* external storage pointer */
  int *OutP;

  /* protect R objects from garbage collector */
  PROTECT(p = AS_NUMERIC(p));
  PROTECT(n = AS_INTEGER(n));
  PROTECT(guess = AS_INTEGER(guess));

  /* Translate arguments from R to C */
  N = *(INTEGER_POINTER(n));
  P = *(NUMERIC_POINTER(p));
  nwmax = *(INTEGER_POINTER(guess));

  /* Allocate space for result */
  w = (int *) R_alloc(nwmax, sizeof(int));

  /* set up */
  GetRNGstate();
  log1p = -log(1.0 - P);
  
  /* main loop */
  i = 0;  /* last selected element of 1...N */
  nw = 0;  /* number of selected elements */
  while(i <= N) {
    log1u = exp_rand();  /* an exponential rv is equivalent to -log(1-U) */
    j = (int) ceil(log1u/log1p); /* j is geometric(p) */
    i += j;
    if(nw >= nwmax) {
      /* overflow; allocate more space */
      w  = (int *) S_realloc((char *) w,  2 * nwmax, nwmax, sizeof(int));
      nwmax    = 2 * nwmax;
    }
    /* add 'i' to output vector */
    w[nw] = i;
    ++nw;
  }
  /* The last saved 'i' could have exceeded 'N' */
  /* For efficiency we don't check this in the loop */
  if(nw > 0 && w[nw-1] > N) 
    --nw;

  PutRNGstate();

  /* create result vector */
  PROTECT(Out = NEW_INTEGER(nw));

  /* copy results into output */
  OutP  = INTEGER_POINTER(Out);
  for(k = 0; k < nw; k++)
    OutP[k] = w[k];

  UNPROTECT(4);
  return(Out);
}
