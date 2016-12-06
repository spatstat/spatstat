/*

  linknnd.h

  k-th nearest neighbours in a linear network

  Using sparse representation of network
  ! Data points must be ordered by segment index !

  This code is #included several times in linknnd.c

  Macros required:
  FNAME   Function name
  CROSS   #defined for X-to-Y, undefined for X-to-X
  HUH     debugging flag

  $Revision: 1.2 $  $Date: 2016/12/04 12:34:19 $

 */

#define MAT(MATRIXNAME, INDEX, ORDER) MATRIXNAME[(ORDER) + (INDEX) * Kmax]
#define NNDIST(INDEX, ORDER) MAT(nndist, (INDEX), (ORDER))
#define NNWHICH(INDEX, ORDER) MAT(nnwhich, (INDEX), (ORDER))
#define VDIST(INDEX, ORDER) MAT(dminvert, (INDEX), (ORDER))
#define VWHICH(INDEX, ORDER) MAT(whichvert, (INDEX), (ORDER))

#define UPDATENN(INDEX, D, J)	\
  UpdateKnnList(D, J, \
		nndist + (INDEX) * Kmax, \
		nnwhich + (INDEX) * Kmax, \
		Kmax, \
		(double) 0.0)

/* ................. */

void FNAME(kmax,         /* number of neighbours required */
	   np, sp, tp,   /* source data points (ordered by sp) */
#ifdef CROSS
	   nq, sq, tq,   /* target data points (ordered by sq) */
#endif
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
  int *np, *nv, *ns;  /* number of points, vertices, segments */
  int *sp, *from, *to; /* integer vectors (mappings) */
  double *tp; /* fractional location coordinates */
#ifdef CROSS
  int *nq, *sq;
  double *tq;
#endif
  double *huge, *tol;
  double *seglen;
  double *nndist;
  int *nnwhich;
{
  int Np, Nv, Kmax, Nout, i, j, ivleft, ivright, jfirst, jlast, k, m;
  double d, hugevalue, slen, tpi, deltad;
  double *dminvert;  /* min dist from each vertex */
  int *whichvert;   /* which min from each vertex */
  int linvknndist(), UpdateKnnList();

#ifdef CROSS
  int Nq;
#else 
#define Nq Np
#define nq np
#define sq sp
#define tq tp
#endif

  Kmax = *kmax;
  Np = *np;
  Nv = *nv;
  hugevalue = *huge;

#ifdef CROSS
  Nq = *nq;
#endif

  /* First compute min distances to target set from each vertex */
#ifdef HUH
  Rprintf("Computing distances from each vertex\n");
#endif

  dminvert = (double *) R_alloc(Nv * Kmax, sizeof(double));
  whichvert = (int *) R_alloc(Nv * Kmax, sizeof(int));

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

#undef MAT
#undef NNDIST
#undef NNWHICH
#undef VDIST
#undef VWHICH
#undef UPDATENN

#ifndef CROSS
#undef nq
#undef Nq
#undef sq
#undef tq
#endif
