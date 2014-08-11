/*
       distmapbin.c

       Distance transform of a discrete binary image
       (8-connected path metric)
       
       $Revision: 1.6 $ $Date: 2011/11/20 03:34:16 $

       
*/

#include <math.h>
#include "raster.h"
#include <R_ext/Utils.h>

void   dist_to_bdry();
void   shape_raster();


void
distmap_bin(in, dist)
        Raster  *in;            /* input:  binary image */
	Raster	*dist;		/* output: distance to nearest point */
	/* rasters must have been dimensioned by shape_raster()
	   and must all have identical dimensions and margins */
{
	int	j,k;
	double	d, dnew;
	double  xstep, ystep, diagstep, huge;
	int rmin, rmax, cmin, cmax;

	/* distances between neighbouring pixels */
	xstep = in->xstep;
	ystep = in->ystep;
	diagstep = sqrt(xstep * xstep + ystep * ystep);
	if(xstep < 0) xstep = -xstep;
	if(ystep < 0) ystep = -ystep;

	/* effectively infinite distance */
	huge = 2.0 * Distance(dist->xmin,dist->ymin,dist->xmax,dist->ymax); 

	/* image boundaries */
	rmin = in->rmin;
	rmax = in->rmax;
	cmin = in->cmin;
	cmax = in->cmax;

#define DISTANCE(ROW, COL) Entry(*dist, ROW, COL, double)
#define MASKTRUE(ROW, COL) (Entry(*in, ROW, COL, int) != 0)
#define MASKFALSE(ROW, COL) (Entry(*in, ROW, COL, int) == 0)
#define UPDATE(D, ROW, COL, STEP) \
	dnew = STEP + DISTANCE(ROW, COL); \
        if(D > dnew) D = dnew

	/* initialise edges to boundary condition */
	for(j = rmin-1; j <= rmax+1; j++) {
	  DISTANCE(j, cmin-1) = (MASKTRUE(j, cmin-1)) ? 0.0 : huge;
	  DISTANCE(j, cmax+1) = (MASKTRUE(j, cmax+1)) ? 0.0 : huge;
	}
	for(k = cmin-1; k <= cmax+1; k++) {
	  DISTANCE(rmin-1, k) = (MASKTRUE(rmin-1, k)) ? 0.0 : huge;
	  DISTANCE(rmax+1, k) = (MASKTRUE(rmax+1, k)) ? 0.0 : huge;
	}
	  
	/* forward pass */

	for(j = rmin; j <= rmax; j++) {
	  R_CheckUserInterrupt();
	  for(k = cmin; k <= cmax; k++) {
	    if(MASKTRUE(j, k))
	      d = DISTANCE(j, k) = 0.0;
	    else {
	      d = huge;
	      UPDATE(d, j-1, k-1, diagstep);
	      UPDATE(d, j-1,   k, ystep);
	      UPDATE(d, j-1, k+1, diagstep);
	      UPDATE(d,   j, k-1, xstep);
	      DISTANCE(j,k) = d;
	    }
	  }
	}

	/* backward pass */

	for(j = rmax; j >= rmin; j--) {
	  R_CheckUserInterrupt();
	  for(k = cmax; k >= cmin; k--) {
	    if(MASKFALSE(j,k)) {
	      d = DISTANCE(j,k);
	      UPDATE(d, j+1, k+1, diagstep);
	      UPDATE(d, j+1,   k, ystep);
	      UPDATE(d, j+1, k-1, diagstep);
	      UPDATE(d,   j, k+1, xstep);
	      DISTANCE(j,k) = d;
	    } 
	  }
	}
}

/* R interface */

void distmapbin(xmin, ymin, xmax, ymax, nr, nc,
	   in, distances, boundary)
	double *xmin, *ymin, *xmax, *ymax;  	  /* x, y dimensions */
	int *nr, *nc;	 	                  /* raster dimensions
				                     EXCLUDING margin of 1 on each side */
	int   *in;              /* input:  binary image */
	double *distances;	/* output: distance to nearest point */
	double *boundary;       /* output: distance to boundary of rectangle */
	/* all images must have identical dimensions including a margin of 1 on each side */
{
	Raster data, dist, bdist;

	shape_raster( &data, (void *) in, *xmin,*ymin,*xmax,*ymax,
			    *nr+2, *nc+2, 1, 1);
	shape_raster( &dist, (void *) distances,*xmin,*ymin,*xmax,*ymax,
			   *nr+2,*nc+2,1,1);
	shape_raster( &bdist, (void *) boundary, *xmin,*ymin,*xmax,*ymax,
			   *nr+2,*nc+2,1,1);
	
	distmap_bin(&data, &dist);

	dist_to_bdry(&bdist);
}	
