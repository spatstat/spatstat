/*
       exactPdist.c

       `Pseudoexact' distance transform of a discrete binary image
       (the closest counterpart to `exactdist.c')
       
       $Revision: 1.13 $ $Date: 2018/12/18 02:43:11 $

       
  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#include <math.h>
#include "raster.h"

void   dist_to_bdry();
void   shape_raster();

void
ps_exact_dt(in, dist, row, col)
        Raster  *in;            /* input:  binary image */
	Raster	*dist;		/* output: exact distance to nearest point */
	Raster	*row;		/* output: row index of closest point */
	Raster	*col;		/* output: column index of closest point */
	/* rasters must have been dimensioned by shape_raster()
	   and must all have identical dimensions and margins */
{
	int	j,k;
	double	d, x, y;
	int	r, c;
	double	dnew;
	double  huge;
	/*	double  bdiag; */
	
	    /* initialise */
#define UNDEFINED -1
#define Is_Defined(I) (I >= 0)
#define Is_Undefined(I) (I < 0)
	
	Clear(*row,int,UNDEFINED)
	Clear(*col,int,UNDEFINED)
		
	huge = 2.0 * DistanceSquared(dist->xmin,dist->ymin,dist->xmax,dist->ymax); 
	Clear(*dist,double,huge)


	  /* if input pixel is TRUE, set distance to 0 and make pixel point to itself */
	for(j = in->rmin; j <= in->rmax; j++)
	for(k = in->cmin; k <= in->cmax; k++) 
	  if(Entry(*in, j, k, int) != 0) {
	      Entry(*dist, j, k, double) = 0.0;
	      Entry(*row,  j, k, int)   = j;
	      Entry(*col,  j, k, int)   = k;
	  }

	/* how to update the distance values */
	
#define GETVALUES(ROW,COL) \
	x = Xpos(*in, COL); \
	y = Ypos(*in, ROW); \
	d = Entry(*dist,ROW,COL,double); 

#define COMPARE(ROW,COL,RR,CC) \
	r = Entry(*row,RR,CC,int); \
	c = Entry(*col,RR,CC,int); \
	if(Is_Defined(r) && Is_Defined(c) \
	   && Entry(*dist,RR,CC,double) < d) { \
	     dnew = DistanceSquared(x, y, Xpos(*in,c), Ypos(*in,r)); \
	     if(dnew < d) { \
		Entry(*row,ROW,COL,int) = r; \
		Entry(*col,ROW,COL,int) = c; \
		Entry(*dist,ROW,COL,double) = dnew; \
		d = dnew; \
	     } \
	}

	/* bound on diagonal step distance squared */
	/* bdiag = (in->xstep * in->xstep + in->ystep * in->ystep); */
	
	/* forward pass */

	for(j = in->rmin; j <= in->rmax; j++)
	for(k = in->cmin; k <= in->cmax; k++) {
	        GETVALUES(j, k)
		COMPARE(j,k, j-1,k-1)
		COMPARE(j,k, j-1,  k)
		COMPARE(j,k, j-1,k+1)
		COMPARE(j,k, j,  k-1)
		  }

	/* backward pass */

	for(j = in->rmax; j >= in->rmin; j--) 
	for(k = in->cmax; k >= in->cmin; k--) {
	        GETVALUES(j, k)
		COMPARE(j,k, j+1,k+1)
		COMPARE(j,k, j+1,  k)
		COMPARE(j,k, j+1,k-1)
		COMPARE(j,k, j,  k+1)
		  } 

	/* take square roots of distances^2 */

	for(j = in->rmax; j >= in->rmin; j--) 
	for(k = in->cmax; k >= in->cmin; k--) 
	        Entry(*dist,j,k,double) = sqrt(Entry(*dist,j,k,double));

}

/* R interface */

void ps_exact_dt_R(xmin, ymin, xmax, ymax, nr, nc, mr, mc, 
	   inp, distances, rows, cols, boundary)
	double *xmin, *ymin, *xmax, *ymax;  	  /* x, y dimensions */
	int *nr, *nc;	 	                  /* raster dimensions
				                     EXCLUDING margins */
	int *mr, *mc;                             /* margins */
	int   *inp;              /* input:  binary image */
	double *distances;	/* output: distance to nearest point */
	int   *rows;	        /* output: row of nearest point (start= 0) */
	int   *cols;	        /* output: column of nearest point (start = 0) */
	double *boundary;       /* output: distance to boundary of rectangle */
	/* all images must have identical dimensions including a margin of 1 on each side */
{
	Raster data, dist, row, col, bdist;
	int mrow, mcol, nrow, ncol;

	mrow = *mr;
	mcol = *mc;

	/* full dimensions */
	nrow = *nr + 2 * mrow;
	ncol = *nc + 2 * mcol;

	shape_raster( &data, (void *) inp, *xmin,*ymin,*xmax,*ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &dist, (void *) distances, *xmin,*ymin,*xmax,*ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &row, (void *) rows, *xmin,*ymin,*xmax,*ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &col, (void *) cols, *xmin,*ymin,*xmax,*ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &bdist, (void *) boundary, *xmin,*ymin,*xmax,*ymax,
		      nrow, ncol, mrow, mcol);
	
	ps_exact_dt(&data, &dist, &row, &col);

	dist_to_bdry(&bdist);
}	
