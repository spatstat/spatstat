/* 
   linnncross.h

   Function body definitions with macros

   $Revision: 1.2 $  $Date: 2015/11/28 02:02:50 $

   Macros used:
   FNAME   name of function
   EXCLU   whether serial numbers are provided
   WHICH   whether 'nnwhich' is required

*/

void 
FNAME(np, xp, yp,   /* data points 'from' */
      nq, xq, yq,   /* data points 'to' */
      nv, xv, yv,   /* network vertices */
      ns, from, to,  /* segments */
      dpath,  /* shortest path distances between vertices */
      psegmap, /* map from data points to segments */
      qsegmap, /* map from data points to segments */
#ifdef EXCLU
      idP, idQ, /* serial numbers for patterns p and q */
#endif
      huge, /* value taken as infinity */
      /* OUTPUT */
#ifdef WHICH
      nndist,  /* nearest neighbour distance for each point */
      nnwhich  /* identifies nearest neighbour */
#else 
      nndist  /* nearest neighbour distance for each point */
#endif
)
  int *np, *nq, *nv, *ns;
int *from, *to, *psegmap, *qsegmap; /* integer vectors (mappings) */
#ifdef EXCLU
  int *idP, *idQ;
#endif
  double *xp, *yp, *xq, *yq, *xv, *yv; /* vectors of coordinates */
  double *huge;
  double *dpath; /* matrix */
  double *nndist; /* nearest neighbour distance for each point */
#ifdef WHICH
  int *nnwhich; /* identifies nearest neighbour */
#endif
{
  int Np, Nq, Nv, i, j;
  int segPi, segQj, nbi1, nbi2, nbj1, nbj2; 
  double d, xpi, ypi, xqj, yqj, dXi1, dXi2, d1Xj, d2Xj, d11, d12, d21, d22; 
  double dmin, hugevalue;
#ifdef EXCLU
  int idPi;
#endif
#ifdef WHICH
  int whichmin;
#endif 

  Np = *np;
  Nq = *nq;
  Nv = *nv;
  hugevalue = *huge;

  /* initialise nn distances */
  for(i = 0; i < Np; i++) {
    nndist[i] = hugevalue;
#ifdef WHICH
    nnwhich[i] = -1;
#endif
  }

  /* main loop */
  for(i = 0; i < Np; i++) {
    xpi = xp[i];
    ypi = yp[i];
#ifdef EXCLU
    idPi = idP[i];
#endif
    segPi = psegmap[i];
    nbi1 = from[segPi];
    nbi2 = to[segPi];
    dXi1 = EUCLID(xpi, ypi, xv[nbi1], yv[nbi1]);
    dXi2 = EUCLID(xpi, ypi, xv[nbi2], yv[nbi2]);
    dmin = nndist[i];
#ifdef WHICH
    whichmin = nnwhich[i];
#endif
    for(j = 0; j < Nq; j++) {
#ifdef EXCLU
      if(idQ[j] != idPi) {
#endif
      xqj = xq[j];
      yqj = yq[j];
      segQj = qsegmap[j];
      /* compute path distance between i and j */
      if(segPi == segQj) {
	/* points i and j lie on the same segment; use Euclidean distance */
	d = EUCLID(xpi, ypi, xqj, yqj);
      } else {
	/* Shortest path from i to j passes through ends of segments;
	   Calculate shortest of 4 possible paths from i to j
	*/
	nbj1 = from[segQj];
	nbj2 = to[segQj];
	d1Xj = EUCLID(xv[nbj1], yv[nbj1], xqj, yqj);
	d2Xj = EUCLID(xv[nbj2], yv[nbj2], xqj, yqj);
	d11 = dXi1 + DPATH(nbi1,nbj1) + d1Xj;
	d12 = dXi1 + DPATH(nbi1,nbj2) + d2Xj;
	d21 = dXi2 + DPATH(nbi2,nbj1) + d1Xj;
	d22 = dXi2 + DPATH(nbi2,nbj2) + d2Xj;
	d = d11;
	if(d12 < d) d = d12;
	if(d21 < d) d = d21;
	if(d22 < d) d = d22;
      }
      /* OK, distance between i and j is d */

      /* update nn for point i */
      if(d < dmin) {
	dmin = d;
#ifdef WHICH
	whichmin = j;
#endif
      }
#ifdef EXCLU
    }
#endif
    }
    /* commit nn distance for point i */
    nndist[i] = dmin;
#ifdef WHICH
    nnwhich[i] = whichmin;
#endif
  }
}

