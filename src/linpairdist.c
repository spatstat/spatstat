#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/* 
   linpairdist.c

   Shortest-path distances between each pair of points in linear network

   $Revision: 1.5 $  $Date: 2012/10/12 10:21:46 $

   linpairdist

*/

#define DPATH(I,J) dpath[(I) + Nv * (J)]
#define ANSWER(I,J) answer[(I) + Np * (J)]
#define EUCLID(X,Y,U,V) sqrt(pow((X)-(U),2)+pow((Y)-(V),2))

void 
linpairdist(np, xp, yp,   /* data points */
	    nv, xv, yv,   /* network vertices */
	    ns, from, to,  /* segments */
	    dpath,  /* shortest path distances between vertices */
	    segmap, /* map from data points to segments */
	    /* OUTPUT */
	    answer  /* shortest path distances between points */
)
  int *np, *nv, *ns;
  int *from, *to, *segmap; /* integer vectors (mappings) */
  double *xp, *yp, *xv, *yv; /* vectors of coordinates */
  double *dpath, *answer; /* matrices */
{
  int Np, Nv, i, j, Np1, maxchunk;
  int segi, segj, nbi1, nbi2, nbj1, nbj2; 
  double d, xpi, ypi, xpj, ypj, dXi1, dXi2, d1Xj, d2Xj, d11, d12, d21, d22; 

  Np = *np;
  Nv = *nv;
  Np1 = Np - 1;

  OUTERCHUNKLOOP(i, Np1, maxchunk, 1024) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Np1, maxchunk, 1024) {
      xpi = xp[i];
      ypi = yp[i];
      segi = segmap[i];
      nbi1 = from[segi];
      nbi2 = to[segi];
      dXi1 = EUCLID(xpi, ypi, xv[nbi1], yv[nbi1]);
      dXi2 = EUCLID(xpi, ypi, xv[nbi2], yv[nbi2]);
      for(j = i+1; j < Np; j++) {
	xpj = xp[j];
	ypj = yp[j];
	segj = segmap[j];
	if(segi == segj) {
	  /* points i and j lie on the same segment; use Euclidean distance */
	  d = sqrt(pow(xpi - xpj, 2) + pow(ypi - ypj, 2));
	} else {
	  /* Shortest path from i to j passes through ends of segments;
	     Calculate shortest of 4 possible paths from i to j
	  */
	  nbj1 = from[segj];
	  nbj2 = to[segj];
	  d1Xj = EUCLID(xv[nbj1], yv[nbj1], xpj, ypj);
	  d2Xj = EUCLID(xv[nbj2], yv[nbj2], xpj, ypj);
	  d11 = dXi1 + DPATH(nbi1,nbj1) + d1Xj;
	  d12 = dXi1 + DPATH(nbi1,nbj2) + d2Xj;
	  d21 = dXi2 + DPATH(nbi2,nbj1) + d1Xj;
	  d22 = dXi2 + DPATH(nbi2,nbj2) + d2Xj;
	  d = d11;
	  if(d12 < d) d = d12;
	  if(d21 < d) d = d21;
	  if(d22 < d) d = d22;
	}
	/* write */
	ANSWER(i,j) = ANSWER(j,i) = d;
      }
      ANSWER(i,i) = 0;
    }
  }
}

