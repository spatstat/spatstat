#include <R.h>
#include "yesno.h"

/* 
   linScrossdist.c

   Distances between points on a linear network
   One pattern to another pattern

   linScrossdist

   'Sparse version' 

   $Revision: 1.4 $  $Date: 2020/04/27 00:52:04 $

   Works with sparse representation
   Requires point data to be ordered by segment index.

   Macros used:
   VERBOSE    debugging

   ! Data points must be ordered by segment index !

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#undef VERBOSE

void Clinvdist();  /* function from linvdist.c */

void 
linScrossdist(np, sp, tp,  /* data points 'from' (ordered by sp) */
	      nq, sq, tq, /* data points 'to'   (ordered by sq) */
	      nv, /* number of network vertices */
	      ns, from, to,  /* segments */
	      seglen,  /* segment lengths */
	      huge, /* value taken as infinity */
	      tol, /* tolerance for updating distances */
	      /* OUTPUT */
	      dist  /* matrix of distances from i to j */
)
  int *np, *nq, *nv, *ns;
  int *from, *to, *sp, *sq; /* integer vectors (mappings) */
  double *tp, *tq; /* fractional location coordinates */
  double *huge, *tol;
  double *seglen; 
  double *dist; 
{
  int Np, Nq, Nv, Npairs, i, j, ivleft, ivright, spi, sqj;
  double dleft, dright, dij, hugevalue, slen, tpi, tqj;
  double *dminvert;  /* min dist from each vertex */
  int one;

  Np = *np;
  Nq = *nq;
  Nv = *nv;
  hugevalue = *huge;
  Npairs = Np * Nq;

  one = 1;

  dminvert = (double *) R_alloc(Nv, sizeof(double));

#ifdef VERBOSE
  Rprintf("Start loop through target points j\n");
#endif

  for(j = 0; j < Nq; j++) {
    R_CheckUserInterrupt();
      
    sqj = sq[j];
    tqj = tq[j];
      
#ifdef VERBOSE
    Rprintf("Target point %d\n\t lies on segment %d\n", j, sqj);
    Rprintf("Compute distance to target from each vertex..\n");
#endif
      
    /* First compute min distance to target point j from each vertex */
    Clinvdist(&one, sq+j, tq+j,
	      nv, ns, from, to, seglen, huge, tol,
	      dminvert);

#ifdef VERBOSE
    Rprintf("Run through source points..\n");
#endif
  
    for(i = 0; i < Np; i++) {
      tpi = tp[i];
      spi = sp[i];   /* segment containing this point */
      slen = seglen[spi];
      
      if(spi == sqj) {
	/* target point lies in the same segment */
	dij = slen * fabs(tqj - tpi);
#ifdef VERBOSE
	Rprintf("\tSource %d and target lie on same segment, distance %lf\n",
		i, dij);
#endif
      } else {
	ivleft = from[spi];
	ivright = to[spi];
#ifdef VERBOSE
	Rprintf("\tSource point %d lies on segment %d = [%d,%d]\n", 
		i, spi, ivleft, ivright);
#endif
	dleft  = slen * tpi + dminvert[ivleft];
	dright = slen * (1.0 - tpi) + dminvert[ivright];
	dij = (dleft < dright) ? dleft : dright;
#ifdef VERBOSE
	Rprintf("\tDistance to left endpoint %d is %lf\n", ivleft, dleft);
	Rprintf("\tDistance to right endpoint %d is %lf\n", ivright, dright);
#endif
      }
#ifdef VERBOSE
      Rprintf("\tAssigning distance d[%d, %d] = %lf\n", i, j, dij);
#endif
      dist[i + j * Np] = dij;
    }
  }
}


