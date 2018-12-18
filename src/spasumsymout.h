/*
  spasumsymout.h

  Function definitions for 'sumsymouter' for sparse matrices/arrays

  This file is #included in sparselinalg.c several times.

  Macros used 

  FNAME     function name

  DBG       (#ifdef) debug 

  WEIGHTS   (#ifdef) use weights 

  $Revision: 1.6 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

void FNAME(m, n, 
	   lenx, ix, jx, kx, x, 
	   flip, 
#ifdef WEIGHTS
	   lenw, jw, kw, w,
#endif
	   y) 
  int *m, *n;            /* dimensions of array m * n * n */
  int *lenx;             /* number of nonzero entries in sparse array x */
  int *ix, *jx, *kx;     /* indices of entries in sparse array x */
  double *x;             /* values in sparse array x */
                         /* NB: ix, jx, kx are assumed to be
			    sorted by order(j,k,i)
			    i.e. in increasing order of j, 
			    then k within j, 
			    then i within (j,k) */
  int *flip;             /* reordering of ix, jx, kx, x that would achieve
			    increasing order(k,j,i) */
#ifdef WEIGHTS
  int *lenw;             /* length of jw, kw */
  int *jw, *kw;          /* indices of entries in sparse matrix w of weights */
                         /* Assumed sorted by order (j,k) */
  double *w;             /* values of weights w */
#endif
  double *y;             /* output: full m * m matrix */
{
  /* Compute the sum of outer(x[,j,k], x[,k,j]) for all j != k */
  int M,N,L, i,j,k,ii, l, ll, lstart, lend, t, tstart, tend, r;
  double xijk, xx;
  int *it, *jt, *kt;
  double *xt;
#ifdef WEIGHTS
  int R;
  double wjk;
#endif

  M = *m; 
  N = *n;
  L = *lenx;
#ifdef WEIGHTS
  R = *lenw;
#endif

  if(L <= 1 || N <= 1 || M <= 0) return;

  /* Create space to store array in k-major order*/
  it = (int *) R_alloc(L, sizeof(int));
  jt = (int *) R_alloc(L, sizeof(int));
  kt = (int *) R_alloc(L, sizeof(int));
  xt = (double *) R_alloc(L, sizeof(double));
  /* copy reordered array */
#ifdef DBG
  Rprintf("----------  Reordered: -------------------\n");
#endif
  for(l = 0; l < L; l++) {
    ll = flip[l];
    it[l] = ix[ll];
    jt[l] = jx[ll];
    kt[l] = kx[ll];
    xt[l] = x[ll];
#ifdef DBG
    Rprintf("%d \t [%d, %d, %d] = %lf\n", l, it[l], jt[l], kt[l], xt[l]);
#endif
  }

  /* Now process array */
  lstart = tstart = r = 0;

  lend = tend = -1; /* to keep compiler happy */

  while(lstart < L && tstart < L) {
    /* Consider a new entry x[,j,k] */
    j = jx[lstart];
    k = kx[lstart];
#ifdef DBG
    Rprintf("Entry %d: [, %d, %d]\n", lstart, j, k);
#endif
#ifdef WEIGHTS     
    /* Find weight w[j,k] */
    while(r < R && ((jw[r] < j) || 
		    ((jw[r] == j) && (kw[r] < k))))
      ++r;
    if(r < R && jw[r] == j && kw[r] == k) {
      /* weight w[j,k] is present */
      wjk = w[r];
#endif
      /* Find all entries in x with the same j,k */
      for(lend = lstart+1;
	  lend < L && jx[lend] == j && kx[lend] == k;
	  ++lend) 
	;
      --lend;
#ifdef DBG
      Rprintf("\t lstart=%d, lend=%d\n", lstart, lend);
#endif
      /* Find corresponding entries in transpose (k'=j, j'=k) */
      /* search forward to find start of run */
      while(tstart < L && ((kt[tstart] < j) ||
			   (kt[tstart] == j && jt[tstart] < k)))
	++tstart;
#ifdef DBG
      Rprintf("\t tstart=%d\n", tstart);
      Rprintf("\t kt[tstart]=%d, jt[tstart]=%d\n", kt[tstart], jt[tstart]);
#endif
      if(tstart < L && kt[tstart] == j && jt[tstart] == k) {
	/* Both x[,j,k] and x[,k,j] are present so a contribution will occur */
	/* seek end of run */
	for(tend = tstart+1;
	    tend < L && kt[tend] == j && jt[tend] == k;
	    ++tend) 
	  ;
	--tend;
#ifdef DBG
	Rprintf("\t tend=%d\n", tend);
#endif
	/* Form products */
	for(l = lstart; l <= lend; l++) {
	  i = ix[l];
	  xijk =  x[l];
#ifdef DBG
	  Rprintf("Entry %d: [%d, %d, %d] = %lf\n", l, i, j, k, xijk);
#endif
	  for(t = tstart; t <= tend; t++) {
	    ii = it[t];
	    xx = xijk * xt[t];
#ifdef WEIGHTS
	    xx *= wjk;
#endif
	    /* increment result at [i, ii] and [ii, i]*/
	    y[i + M * ii] += xx;
	    /*	  y[ii + M * i] += xx;  */
#ifdef DBG
	    Rprintf("-- matches entry %d: [%d, %d, %d] = %lf\n", 
		    t, ii, k, j, xt[t]);
	    Rprintf("++ %lf\n", xx);
#endif
	  }
	}
      }
#ifdef WEIGHTS
    }
#endif
    lstart = ((lend > lstart) ? lend : lstart) + 1;
    tstart = ((tend > tstart) ? tend : tstart) + 1;
  }
}

