#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/* 
   linearradius.c

   Bounding radius in linear network

   $Revision: 1.2 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#define DPATH(I,J) dpath[(J) + Nv * (I)]

#include "yesno.h"

#undef DEBUG

void 
linearradius(ns, from, to,  /* network segments */
	     lengths, /* segment lengths */
	     nv, dpath,  /* shortest path distances between vertices */
	     huge, 
	     result)
     int *nv, *ns;
     int *from, *to; /* integer vectors (mappings) */
     double *dpath; /* matrix of shortest path distances between vertices */
     double *lengths; /* vector of segment lengths */
     double *huge; /* very large value */
     double *result; 
{
  int Nv, Ns;
  int i, j, A, B, C, D;
  double AB, AC, AD, BC, BD, CD;
  double sAij, sBij, sAiMax, sBiMax, smin;
  int maxchunk;

  Nv = *nv;
  Ns = *ns;
  smin = *huge;

  OUTERCHUNKLOOP(i, Ns, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Ns, maxchunk, 16384) {
      /* indices of endpoints of segment i */
      A = from[i];
      B = to[i];
      AB = lengths[i];
      sAiMax = sBiMax = AB/2.0;
      for(j = 0; j < Ns; j++) {
	if(j != i) {
	  /* indices of endpoints of segment i */
	  C = from[j];
	  D = to[j];
	  CD = lengths[j];
	  AC = DPATH(A,C);
	  AD = DPATH(A,D);
	  BC = DPATH(B,C);
	  BD = DPATH(B,D);
	  /* max dist from A to any point in segment j */
	  sAij = (AD > AC + CD) ? AC + CD :
 	          (AC > AD + CD) ? AD + CD : (AC + AD + CD)/2.0;
	  /* max dist from B to any point in segment j */
	  sBij = (BD > BC + CD) ? BC + CD : 
  	          (BC > BD + CD) ? BD + CD : (BC + BD + CD)/2.0;
	  /* row-wise maximum */
	  if(sAij > sAiMax) sAiMax = sAij;
	  if(sBij > sBiMax) sBiMax = sBij;
	}
      }
      if(sAiMax < smin) smin = sAiMax;
      if(sBiMax < smin) smin = sBiMax;
    }
  }

  *result = smin;
}

