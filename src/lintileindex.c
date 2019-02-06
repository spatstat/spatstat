#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/*
  lintileindex.c

  Given a tessellation on a linear network,
  compute tile index for each query point

  NOTE: data are assumed to be sorted by segment

  $Revision: 1.2 $ $Date: 2019/02/06 08:36:02 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

void
lintileindex(n, seg, tp,  /* query points, sorted by segment */
	     dfn, dfseg, dft0, dft1, dftile, /* tessellation data, sorted */
	     answer)
  /* query points */
  int *n;               /* number of query points */
  int *seg;             /* which segment contains this query point*/
  double *tp;           /* position along segment */
  /* tessellation pieces */
  int *dfn;             /* number of pieces */
  int *dfseg;           /* which segment contains this piece */
  double *dft0, *dft1;  /* positions of endpoints of this piece */
  int *dftile;          /* which tile the piece belongs to */
  /* output */
  int *answer;          /* which tile the query point belongs to */
{
  int N, M, i, start, finish, j, segi, currentseg, maxchunk;
  double tpi;
  
  N = *n;
  M = *dfn;
  currentseg  = -1;
  start = finish = 0;

  /* 
     answer[] is implicitly initialised to zero,
     which will serve as NA
 */
  
  OUTERCHUNKLOOP(i, N, maxchunk, 1024) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, N, maxchunk, 1024) {
      segi = seg[i];
      tpi  = tp[i];
      if(segi > currentseg) {
        /* advance the bookmark until reaching data for this segment */
	while(start < M && dfseg[start] < segi) ++start;
	if(start == M) {
	  /* Reached end of data for tessellation */
	  /* All remaining results are NA */
	  return;
	}
	currentseg = dfseg[start];
	for(finish = start;
	    finish < M && dfseg[finish] == currentseg;
	    ++finish) ;
	if(finish == M || dfseg[finish] > currentseg) --finish;
      }
      if(currentseg == segi) {
	for(j = start; j <= finish; ++j) {
	  if(dft0[j] <= tpi && tpi <= dft1[j]) {
	    answer[i] = dftile[j];
	    break;
	  }
	}
      }    
    }
  } 
}
