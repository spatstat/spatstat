/*
  whist.c

  Weighted histogram

  Designed for very fine bins

  Cwhist(indices, weights, nbins)

  indices point to bins (range: 0 to nbins-1)
  
  $Revision: 1.4 $  $Date: 2014/08/07 12:52:54 $

*/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

SEXP Cwhist(SEXP indices, SEXP weights, SEXP nbins) {
  int i, j, a, N, M; 
  int *x;
  double *w, *y;
  SEXP result;

  /* =================== Protect R objects from garbage collector ======= */
  PROTECT(indices = AS_INTEGER(indices));
  PROTECT(weights = AS_NUMERIC(weights));
  PROTECT(nbins   = AS_INTEGER(nbins));

  N = LENGTH(indices);
  M = *(INTEGER_POINTER(nbins));

  x = INTEGER_POINTER(indices);
  w = NUMERIC_POINTER(weights);

  PROTECT(result = NEW_NUMERIC(M));
  y =  NUMERIC_POINTER(result);

  for(j = 0; j < M; j++)
    y[j] = 0.0;

  for(i = 0; i < N; i++) {
    j = x[i];
    if(j != NA_INTEGER && R_FINITE(w[i]) && j >= 0 && j < M)
      y[j] += w[i];
  }
  UNPROTECT(4);
  return(result);
}

