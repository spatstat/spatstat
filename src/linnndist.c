#include <R.h>

/* 
   linnndist.c

   Shortest-path distances between nearest neighbours in linear network

   $Revision: 1.2 $  $Date: 2018/12/18 02:43:11 $

   linnndist
   linnnwhich

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#define DPATH(I,J) dpath[(J) + Nv * (I)]
#define ANSWER(I,J) answer[(J) + Np * (I)]
#define EUCLID(X,Y,U,V) sqrt(pow((X)-(U),2)+pow((Y)-(V),2))

void 
linnndist(np, xp, yp,   /* data points */
	  nv, xv, yv,   /* network vertices */
	  ns, from, to,  /* segments */
	  dpath,  /* shortest path distances between vertices */
	  segmap, /* map from data points to segments */
	  huge, /* value taken as infinity */
	  /* OUTPUT */
	  answer  /* nearest neighbour distance for each point */
)
  int *np, *nv, *ns;
  int *from, *to, *segmap; /* integer vectors (mappings) */
  double *xp, *yp, *xv, *yv; /* vectors of coordinates */
  double *huge;
  double *dpath; /* matrix */
  double *answer; /* vector of output values */
{
  int Np, Nv, i, j, Np1;
  int segi, segj, nbi1, nbi2, nbj1, nbj2; 
  double d, xpi, ypi, xpj, ypj, dXi1, dXi2, d1Xj, d2Xj, d11, d12, d21, d22; 
  double dmin, hugevalue;

  Np = *np;
  Nv = *nv;
  Np1 = Np - 1;
  hugevalue = *huge;

  /* initialise nn distances */
  for(i = 0; i < Np; i++)
    answer[i] = hugevalue;

  /* main loop */
  for(i = 0; i < Np1; i++) {
    xpi = xp[i];
    ypi = yp[i];
    segi = segmap[i];
    nbi1 = from[segi];
    nbi2 = to[segi];
    dXi1 = EUCLID(xpi, ypi, xv[nbi1], yv[nbi1]);
    dXi2 = EUCLID(xpi, ypi, xv[nbi2], yv[nbi2]);
    dmin = answer[i];
    for(j = i+1; j < Np; j++) {
      xpj = xp[j];
      ypj = yp[j];
      segj = segmap[j];
      /* compute path distance between i and j */
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
      /* OK, distance between i and j is d */

      /* update nn distance for point i */
      if(d < dmin) dmin = d;
      /* update nn distance for point j */
      if(d < answer[j]) answer[j] = d;
    }
    /* commit nn distance for point i */
    answer[i] = dmin;
  }
}


void 
linnnwhich(np, xp, yp,   /* data points */
	   nv, xv, yv,   /* network vertices */
	   ns, from, to,  /* segments */
	   dpath,  /* shortest path distances between vertices */
	   segmap, /* map from data points to segments */
	   huge, /* value taken as infinity */
	   /* OUTPUT */
	   nndist,  /* nearest neighbour distance for each point */
	   nnwhich  /* identifies nearest neighbour */
)
  int *np, *nv, *ns;
  int *from, *to, *segmap; /* integer vectors (mappings) */
  double *xp, *yp, *xv, *yv; /* vectors of coordinates */
  double *huge;
  double *dpath; /* matrix */
  double *nndist; /* vector of output values */
  int *nnwhich; /* vector of output values */
{
  int Np, Nv, i, j, Np1;
  int segi, segj, nbi1, nbi2, nbj1, nbj2; 
  double d, xpi, ypi, xpj, ypj, dXi1, dXi2, d1Xj, d2Xj, d11, d12, d21, d22; 
  double dmin, hugevalue;
  int whichmin;

  Np = *np;
  Nv = *nv;
  Np1 = Np - 1;
  hugevalue = *huge;

  /* initialise nn distances and identifiers */
  for(i = 0; i < Np; i++) {
    nndist[i] = hugevalue;
    nnwhich[i] = -1;
  }

  /* main loop */
  for(i = 0; i < Np1; i++) {
    xpi = xp[i];
    ypi = yp[i];
    segi = segmap[i];
    nbi1 = from[segi];
    nbi2 = to[segi];
    dXi1 = EUCLID(xpi, ypi, xv[nbi1], yv[nbi1]);
    dXi2 = EUCLID(xpi, ypi, xv[nbi2], yv[nbi2]);
    dmin = nndist[i];
    whichmin = nnwhich[i];
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
      /* OK, distance between i and j is d */

      /* update nn for point i */
      if(d < dmin) {
	dmin = d;
	whichmin = j;
      }
      /* update nn for point j */
      if(d < nndist[j]) {
	nndist[j] = d;
	nnwhich[j] = i;
      }
    }
    /* commit nn for point i */
    nndist[i] = dmin;
    nnwhich[i] = whichmin;
  }
}
