/*
  sumsymouter.h

  Code template for some functions in linalg.c

  $Revision: 1.3 $ $Date: 2013/04/18 11:55:24 $

  Macros used: FNAME = function name,
               WEIGHTED = #defined for weighted version
*/

void FNAME(
  x, 
#ifdef WEIGHTED
  w, 
#endif
  p, n, y
) 
  double *x;    /* p by n by n array */
#ifdef WEIGHTED
  double *w;    /*      n by n matrix (symmetric) */
#endif
  int *p, *n;
  double *y;    /* output matrix p by p, initialised to zero */
{
  int N, P;
  register int i, j, k, m, ijpos, jipos, maxchunk;
  register double *xij, *xji;
#ifdef WEIGHTED
  register double wij;
#endif
  N = *n; 
  P = *p;
  OUTERCHUNKLOOP(i, N, maxchunk, 256) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, N, maxchunk, 256) {
      /* loop over j != i */
      if(i > 0) {
	for(j = 0; j < i; j++) {
	  /* pointers to [i,j] and [j,i] in N*N matrices */
	  ijpos = i + N * j;
	  jipos = j + N * i;
	  /* pointers to x[, i, j] and x[ , j, i] */
	  xij = x + ijpos * P;
	  xji = x + jipos * P;
	  /* outer product */ 
#ifdef WEIGHTED
	  wij = w[ijpos];
#endif
	  for(k = 0; k < P; k++) {
	    for(m = 0; m < P; m++) {
#ifdef WEIGHTED
	      y[m + k * P] += wij * xij[m] * xji[k];
#else
	      y[m + k * P] += xij[m] * xji[k];
#endif
	    }
	  }
	}
      }
      if(i + 1 < N) {
	for(j = i+1; j < N; j++) {
	  /* pointers to [i,j] and [j,i] in N*N matrices */
	  ijpos = i + N * j;
	  jipos = j + N * i;
	  /* pointers to x[, i, j] and x[ , j, i] */
	  xij = x + ijpos * P;
	  xji = x + jipos * P;
	  /* outer product */ 
#ifdef WEIGHTED
	  wij = w[ijpos];
#endif
	  for(k = 0; k < P; k++) {
	    for(m = 0; m < P; m++) {
#ifdef WEIGHTED
	      y[m + k * P] += wij * xij[m] * xji[k];
#else
	      y[m + k * P] +=       xij[m] * xji[k];
#endif
	    }
	  }
	}
      }
      /* end of loop over j */
    }
  }
}

