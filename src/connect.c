/*
       connect.c

       Connected component transforms

       cocoImage:   connected component transform of a discrete binary image
                   (8-connected topology)

       cocoGraph: connected component labels for a discrete graph
                   specified by a list of edges
       
       $Revision: 1.7 $ $Date: 2013/02/24 03:19:10 $

       
*/

#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

#include "raster.h"
void   shape_raster();

#ifndef TRUE
#define TRUE (0 == 0)
#define FALSE (!TRUE)
#endif


/* workhorse function for cocoImage */

void
comcommer(im)
     Raster  *im;            
     /* raster must have been dimensioned by shape_raster() */
     /* Pixel values assumed to be 0 in background, and 
        distinct nonzero integers in foreground */
{
  int	j,k;
  int rmin, rmax, cmin, cmax;
  int label, curlabel, minlabel;
  int nchanged;

  /* image boundaries */
  rmin = im->rmin;
  rmax = im->rmax;
  cmin = im->cmin;
  cmax = im->cmax;

#define ENTRY(ROW, COL) Entry(*im, ROW, COL, int)

#define UPDATE(ROW,COL,BEST,NEW) \
     NEW = ENTRY(ROW, COL); \
     if(NEW != 0 && NEW < BEST) \
       BEST = NEW

  nchanged = 1;

  while(nchanged >0) {
    nchanged = 0;
    R_CheckUserInterrupt();
    for(j = rmin; j <= rmax; j++) {
      for(k = cmin; k <= cmax; k++) {
	curlabel = ENTRY(j, k);
	if(curlabel != 0) {
	  minlabel = curlabel;
	  UPDATE(j-1, k-1, minlabel, label);
	  UPDATE(j-1, k,   minlabel, label);
	  UPDATE(j-1, k+1, minlabel, label);
	  UPDATE(j,   k-1, minlabel, label);
	  UPDATE(j,   k,   minlabel, label);
	  UPDATE(j,   k+1, minlabel, label);
	  UPDATE(j+1, k-1, minlabel, label);
	  UPDATE(j+1, k,   minlabel, label);
	  UPDATE(j+1, k+1, minlabel, label);
	  if(minlabel < curlabel) {
	    ENTRY(j, k) = minlabel;
	    nchanged++;
	  }
	}
      }
    }
  }
}

void cocoImage(mat, nr, nc)
     int   *mat;        /* input:  binary image */
     int *nr, *nc;      /* raster dimensions
			   EXCLUDING margin of 1 on each side */
{
  Raster im;

  shape_raster( &im, (void *) mat, 
		(double) 1, (double) 1,
		(double) *nc, (double) *nr, 
		*nr+2, *nc+2, 1, 1);
  comcommer(&im);
}	

void cocoGraph(nv, ne, ie, je, label, status)
     /* inputs */
     int *nv;         /* number of graph vertices */
     int *ne;         /* number of edges */
     int *ie, *je;    /* vectors of indices of ends of each edge */ 
     /* output */
     int *label;      /* vector of component labels for each vertex */
                      /* Component label is lowest serial number of
			 any vertex in the connected component */
     int *status;          /* 0 if OK, 1 if overflow */
{
  int Nv, Ne, i, j, k, niter, labi, labj, changed;
  
  Nv = *nv;
  Ne = *ne;

  /* initialise labels */
  for(k = 0; k < Nv; k++)
    label[k] = k;

  for(niter = 0; niter < Nv; niter++) {
    R_CheckUserInterrupt();
    changed = FALSE;
    for(k = 0; k < Ne; k++) {
      i = ie[k];
      j = je[k];
      labi = label[i];
      labj = label[j];
      if(labi < labj) {
	label[j] = labi;
	changed = TRUE;
      } else if(labj < labi) {
	label[i] = labj;
	changed = TRUE;
      } 
    }
    if(!changed) {
      /* algorithm has converged */
      *status = 0;
      return;
    }
  }
  /* error exit */   
  *status = 1;
  return;
}
