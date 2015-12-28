#include <R.h>
#include "yesno.h"

/*

  linknnd.c

  k-th nearest neighbours (source points to target points)
  in a linear network

  Sparse representation of network

  $Revision: 1.1 $  $Date: 2015/12/28 02:46:17 $

  ! Data points must be ordered by segment index !

 */

#undef HUH
#define MAT(MATRIXNAME, INDEX, ORDER) MATRIXNAME[(ORDER) + (INDEX) * Kmax]

void linknncross(kmax,         /* number of neighbours required */
		 np, sp, tp,   /* source data points (ordered by sp) */
		 nq, sq, tq,   /* target data points (ordered by sq) */
		 nv,           /* number of network vertices */
		 ns, from, to, /* segments (pairs of vertices) */
		 seglen,       /* segment lengths */
		 huge,         /* value taken as infinity */
		 tol,          /* tolerance for updating distances */
		 /* OUTPUT */
		 nndist,         /* distance from each source point to
				  the nearest, ..., kth nearest target points */
		 nnwhich         /* identifies which target points */
		 )
  int *kmax;
  int *np, *nq, *nv, *ns;  /* number of points, vertices, segments */
  int *sp, *sq, *from, *to; /* integer vectors (mappings) */
  double *tp, *tq; /* fractional location coordinates */
  double *huge, *tol;
  double *seglen;
  double *nndist;
  int *nnwhich;
{
  int Np, Nq, Nv, Ns, Kmax, Nout, i, j, ivleft, ivright, jfirst, jlast, k, m;
  int segPi, segQj, nbi1, nbi2, nbj1, nbj2; 
  double eps, d, dmin, hugevalue, slen, dleft, dright, tpi, tqj, deltad;
  char converged;
  double *dminvert;  /* min dist from each vertex */
  int whichmin;
  int *whichvert;   /* which min from each vertex */
  int linvknndist(), UpdateKnnList();

  Kmax = *kmax;
  Np = *np;
  Nq = *nq;
  Nv = *nv;
  Ns = *ns;
  hugevalue = *huge;
  eps = *tol;

  /* First compute min distances to target set from each vertex */
#ifdef HUH
  Rprintf("Computing distances from each vertex\n");
#endif
  dminvert = (double *) R_alloc(Nv * Kmax, sizeof(double));
  whichvert = (int *) R_alloc(Nv * Kmax, sizeof(int));
#define VDIST(INDEX, ORDER) MAT(dminvert, (INDEX), (ORDER))
#define VWHICH(INDEX, ORDER) MAT(whichvert, (INDEX), (ORDER))
  linvknndist(kmax, nq, sq, tq, nv, ns, from, to, seglen, huge, tol, 
	     dminvert, whichvert);

#ifdef HUH
  Rprintf("Initialise answer\n");
#endif
  /* initialise nn distances from source points */
  Nout = Np * Kmax;
  for(i = 0; i < Nout; i++) {
    nndist[i] = hugevalue;
    nnwhich[i] = -1;
  }

#define NNDIST(INDEX, ORDER) MAT(nndist, (INDEX), (ORDER))
#define NNWHICH(INDEX, ORDER) MAT(nnwhich, (INDEX), (ORDER))

#define UPDATENN(INDEX, D, J)	\
  UpdateKnnList(D, J, \
		nndist + (INDEX) * Kmax, \
		nnwhich + (INDEX) * Kmax, \
		Kmax, \
		(double) 0.0)


  /* run through all source points */
#ifdef HUH
  Rprintf("Run through source points\n");
#endif
  jfirst = 0;
  for(i = 0; i < Np; i++) {
    tpi = tp[i];
    m = sp[i];   /* segment containing this point */
    slen = seglen[m];
    ivleft = from[m];
    ivright = to[m];
#ifdef HUH
    Rprintf("Source point %d lies on segment %d = [%d,%d]\n", 
	    i, m, ivleft, ivright);
#endif
    deltad = slen * tpi;
#ifdef HUH
    Rprintf("\tComparing to left endpoint %d, distance %lf\n", ivleft, deltad);
#endif
    for(k = 0; k < Kmax; k++)
      UPDATENN(i, deltad + VDIST(ivleft, k), VWHICH(ivleft, k));

    deltad = slen * (1.0 - tpi);
#ifdef HUH
   Rprintf("\tComparing to right endpoint %d, distance %lf\n", ivright, deltad);
#endif
    for(k = 0; k < Kmax; k++)
      UPDATENN(i, deltad + VDIST(ivright, k), VWHICH(ivright, k));

    /* find any target points in this segment */
    while(jfirst < Nq && sq[jfirst] < m) jfirst++;
    jlast = jfirst;
    while(jlast < Nq && sq[jlast] == m) jlast++;
    --jlast;
    /* if there are no such points, then jlast < jfirst */
    if(jfirst <= jlast) {
      for(j = jfirst; j <= jlast; j++) {
	d = slen * fabs(tq[j] - tpi);
	UPDATENN(i, d, j);
      }
    }
  }
}

