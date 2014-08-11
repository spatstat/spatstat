/*

  crossloop.h

  Generic code template for loop for cross-close-pairs operations
  collecting contributions to point x_i
  from all points y_j such that ||x_i - y_j|| <= r

  cpp variables used:

       INITIAL_I        code executed at start of 'i' loop       
       CONTRIBUTE_IJ    code executed to compute contribution from j to i
       COMMIT_I         code executed to save total contribution to i

  C variables used:
       int i, j, n1, n2, maxchunk, jleft;
       double x1i, y1i, xleft, dx, dy, d2, rmax, r2max;
       double *x1, *y1, *x2, *y2;

  $Revision: 1.2 $  $Date: 2014/04/02 07:59:10 $

*/

#ifndef CHUNKLOOP_H
#include "chunkloop.h"
#endif

#define CROSSLOOP(INITIAL_I, CONTRIBUTE_IJ, COMMIT_I)           \
  OUTERCHUNKLOOP(i, n1, maxchunk, 65536) {	     	  	\
    R_CheckUserInterrupt();			     	  	\
    INNERCHUNKLOOP(i, n1, maxchunk, 65536) {	     	  	\
						     	  	\
      x1i = x1[i];					        \
      y1i = y1[i];					        \
                                                     	        \
      INITIAL_I;					  	\
                                                     	        \
      jleft = 0;					  	\
							  	\
      /* 						  	\
	 adjust starting point jleft			  	\
      */						  	\
      xleft = x1i - rmax;				  	\
      while((x2[jleft] < xleft) && (jleft+1 < n2))	  	\
	++jleft;					  	\
							  	\
      /* 						  	\
	 process from j = jleft until dx > rmax		  	\
      */						  	\
      for(j=jleft; j < n2; j++) {			  	\
	dx = x2[j] - x1i;				        \
	if(dx > rmax)					  	\
	  break;					  	\
	dy = y2[j] - y1i;				  	\
	d2 = dx * dx + dy * dy;				  	\
	if(d2 <= r2max) {				  	\
	    /* add this (i, j) pair to output */	  	\
	  CONTRIBUTE_IJ;				  	\
	}						  	\
      }							  	\
      COMMIT_I;						  	\
    }							  	\
  }  
