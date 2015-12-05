/*

  linvdist.h

  Distance function at vertices
  (shortest distance from each vertex to a data point)

  Function body definitions with macros

  Sparse representation of network

  $Revision: 1.3 $  $Date: 2015/12/05 07:26:56 $

  Macros used:
  FNAME   name of function
  WHICH   whether 'nnwhich' is required
  HUH     debugging flag

  ! Data points must be ordered by segment index !

*/

void FNAME(np, sp, tp,   /* target data points (ordered by sp) */
	   nv,           /* number of network vertices */
	   ns, from, to, /* segments */
	   seglen,       /* segment lengths */
	   huge,         /* value taken as infinity */
	   tol,          /* tolerance for updating distances */
	   /* OUTPUT */
#ifdef WHICH
	   dist,         /* distance from each vertex to nearest data point */
	   which         /* identifies nearest data point */
#else 
	   dist          /* distance from each vertex to nearest data point */
#endif	   
) 
  int *np, *nv, *ns;  /* number of points, vertices, segments */
  int *sp, *from, *to; /* integer vectors (mappings) */
  double *tp; /* fractional location coordinates */
  double *huge, *tol;
  double *seglen;
  double *dist;
#ifdef WHICH
  int *which;
#endif
{
  int Np, Nv, Ns, i, j, k, segPj, ivleft, ivright;
  double hugevalue, eps, dleft, dright, slen, d, tpj;
  char converged;

  Np = *np;
  Nv = *nv;
  Ns = *ns;
  hugevalue = *huge;
  eps = *tol;

#ifdef HUH
  Rprintf("Initialise dist\n");
#endif
  /* initialise to huge value */
  for(i = 0; i < Nv; i++) {
    dist[i] = hugevalue;
#ifdef WHICH
    which[i] = -1;
#endif
  }

#ifdef HUH
  Rprintf("Run through target points\n");
#endif
  /* assign correct value to endpoints of segments containing target points */
  for(j = 0; j < Np; j++) {
    segPj = sp[j];
    tpj = tp[j];
    slen = seglen[segPj];
    ivleft = from[segPj];
    d = slen * tpj;
    if(d < dist[ivleft]) {
      dist[ivleft] = d;
#ifdef WHICH
      which[ivleft] = j;
#endif
    }
    ivright = to[segPj];
    d = slen * (1.0 - tpj);
    if(d < dist[ivright]) {
      dist[ivright] = d;
#ifdef WHICH
      which[ivright] = j;
#endif
    }
  }

  /* recursively update */
#ifdef HUH
  Rprintf("Recursive update\n");
#endif
  converged = NO;
  while(!converged) {
    converged = YES;
#ifdef HUH
    Rprintf("........... starting new pass ...................... \n");
#endif
    for(k = 0; k < Ns; k++) {
      ivleft = from[k];
      ivright = to[k];
      slen = seglen[k];
      dleft = (double) dist[ivleft];
      dright = (double) dist[ivright];
      d = (double) (dleft + slen);
      if(d < dright - eps) {
#ifdef HUH
	Rprintf("Updating ivright=%d using ivleft=%d, from %lf to %lf+%lf=%lf\n",
		ivright, ivleft, dright, dleft, slen, d);
#endif
	converged = NO;
	dist[ivright] = d;
#ifdef WHICH
	which[ivright] = which[ivleft];
#endif
      } else {
	d = (double) (dright + slen);
	if(d < dleft - eps) {
#ifdef HUH
	Rprintf("Updating ivleft=%d using ivright=%d, from %lf to %lf+%lf=%lf\n",
		ivleft, ivright, dleft, dright, slen, d);
#endif
	  converged = NO;
	  dist[ivleft] = d;
#ifdef WHICH
	  which[ivleft] = which[ivright];
#endif
	}
      }
    }
  }

#ifdef HUH
  Rprintf("Done\nVertex values:\n");
#ifdef WHICH
  Rprintf("\ti\twhich\tdist\n");
  for(i = 0; i < Nv; i++) 
    Rprintf("\t%d\t%d\t%lf\n", i, which[i], dist[i]);
#else
  Rprintf("\ti\tdist\n");
  for(i = 0; i < Nv; i++) 
    Rprintf("\t%d\t%lf\n", i, dist[i]);
#endif
#endif
}
