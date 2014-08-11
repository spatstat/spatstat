/*

  Function body for dist2dpath.c

  Macros used: 

  FNAME   function name
  DTYPE   declaration for distance values ('double' or 'int')
  FLOATY  (DTYPE == 'double')

  $Revision: 1.2 $   $Date: 2012/10/22 07:31:54 $

 */

#undef DEBUG 

#define MATRIX(X,I,J) (X)[(J) + n * (I)]
#define D(I,J)     MATRIX(d,     I, J)
#define DPATH(I,J) MATRIX(dpath, I, J)
#define ADJ(I,J)   (MATRIX(adj,  I, J) != 0)

#define INFIN -1
#define FINITE(X) ((X) >= 0)

void FNAME(nv, d, adj, dpath, tol, niter, status) 
  int *nv;     /* number of vertices */
  DTYPE *d;  /* matrix of edge lengths */
  int *adj;   /* 0/1 edge matrix of graph */
  DTYPE *tol;  /* tolerance threshold (ignored in integer case) */
  DTYPE *dpath; /* output - shortest path distance matrix */
  int *niter, *status; /* status = 0 for convergence */
{
  int i, j, k, n, iter, maxiter, changed;
  DTYPE dij, dik, dkj, dikj;
#ifdef FLOATY
  DTYPE eps, diff, maxdiff;
#endif
  int totaledges, starti, nneighi, increm, pos;
  int *start, *nneigh, *indx;

  n = *nv;
#ifdef FLOATY
  eps = *tol;
#endif

  /* initialise and count edges */
  *status = -1;
  totaledges = 0;
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      DPATH(i, j) = (i == j) ? 0 : ((ADJ(i,j)) ? D(i, j) : INFIN);
      if((i != j) && ADJ(i,j)) ++totaledges;
    }
  }

  maxiter = 2 + ((totaledges > n) ? totaledges : n);

  /* store indices j for each edge (i,j) */
  indx = (int *) R_alloc(totaledges, sizeof(int));
  nneigh = (int *) R_alloc(n, sizeof(int));
  start  = (int *) R_alloc(n, sizeof(int));

  pos = 0;
  for(i = 0; i < n; i++) {
    nneigh[i] = 0;
    start[i] = pos;
#ifdef DEBUG 
    Rprintf("Neighbours of %d:\n", i);
#endif
    for(j = 0; j < n; j++) {
      if((i != j) && ADJ(i,j) && FINITE(D(i,j))) {
#ifdef DEBUG 
	Rprintf("\t%d\n", j);
#endif
	++(nneigh[i]);
	if(pos > totaledges)
	  error("internal error: pos exceeded storage");
	indx[pos] = j;
	++pos;
      }
    }
  }

  /* run */
  for(iter = 0; iter < maxiter; iter++) {

    changed = 0;
#ifdef FLOATY
    maxdiff = 0;
#endif

#ifdef DEBUG
    Rprintf("--------- iteration %d ---------------\n", iter);
#endif
    for(i = 0; i < n; i++) {
      R_CheckUserInterrupt();
      nneighi = nneigh[i];
      if(nneighi > 0) {
	/* run through neighbours k of i */
	starti = start[i];
	for(increm = 0, pos=starti; increm < nneighi; ++increm, ++pos) {
	  k = indx[pos];
	  dik = DPATH(i,k);
#ifdef DEBUG
#ifdef FLOATY
	    Rprintf("i=%d k=%d dik=%lf\n", i, k, dik);
#else
	    Rprintf("i=%d k=%d dik=%d\n",  i, k, dik);
#endif
#endif
	  /* now run through all other vertices j */
	  for(j = 0; j < n; j++) {
	    if(j != i && j != k) {
	      dij = DPATH(i,j);
	      dkj = DPATH(k,j);
	      if(FINITE(dkj)) {
		dikj = dik + dkj;
#ifdef DEBUG
#ifdef FLOATY
		Rprintf("considering %d -> (%d) -> %d,\t dij=%lf, dikj=%lf\n", 
			i, k, j, dij, dikj);
#else
		Rprintf("considering %d -> (%d) -> %d,\t dij=%d, dikj=%d\n", 
			i, k, j, dij, dikj);
#endif
#endif
		if(!FINITE(dij) || dikj < dij) {
#ifdef DEBUG
#ifdef FLOATY
		  Rprintf("updating i=%d j=%d via k=%d from %lf to %lf\n", 
			  i, j, k, dij, dikj);
#else
		  Rprintf("updating i=%d j=%d via k=%d from %d to %d\n", 
			  i, j, k, dij, dikj);
#endif
#endif
		  DPATH(i,j) = DPATH(j,i) = dikj;
		  changed = 1;
#ifdef FLOATY
		  diff = (FINITE(dij)) ? dij - dikj : dikj;
		  if(diff > maxdiff) maxdiff = diff;
#endif
		}
	      }
	    }
	  }
	}
      }
    }
    if(changed == 0) {
      /* algorithm converged */
#ifdef DEBUG
      Rprintf("Algorithm converged\n");
#endif
      *status = 0;
      break;
#ifdef FLOATY
    } else if(FINITE(maxdiff) && maxdiff < eps) {
      /* tolerance reached */
#ifdef DEBUG
      Rprintf("Algorithm terminated with maxdiff=%lf\n", maxdiff);
#endif
      *status = 1;
      break;
#endif
    }
  }

#ifdef DEBUG
  Rprintf("Returning after %d iterations on %d vertices\n", iter, n);
#endif
  
  *niter = iter;
}

#undef DEBUG 

#undef MATRIX
#undef D
#undef DPATH
#undef ADJ
#undef INFIN 
#undef FINITE

