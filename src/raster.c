/*
  
  raster.c

  shape_raster()     initialise a Raster structure
  
  $Revision: 1.1 $ $Date: 2020/11/30 11:19:18 $

*/

#include <math.h>
#include "raster.h"

void 
shape_raster(ras,data,xmin,ymin,xmax,ymax,nrow,ncol,mrow,mcol)
     Raster          *ras;           /* the raster structure to be initialised */
     void		*data;
     int 	        nrow, ncol;  /* absolute dimensions of storage array */
     int 		mrow, mcol;  /* margins clipped off */
	                             /* e.g. valid width is ncol - 2*mcol columns */
     double		xmin, ymin,	/* image dimensions in R^2 after clipping */
		        xmax, ymax;     
{
	ras->data	= data;
	ras->nrow 	= nrow;
	ras->ncol 	= ncol;
	ras->length 	= nrow * ncol;
	ras->rmin	= mrow;
	ras->rmax	= nrow - mrow - 1;
	ras->cmin	= mcol;
	ras->cmax	= ncol - mcol - 1;
	ras->x0		= 
	ras->xmin	= xmin;
	ras->x1 	=
	ras->xmax	= xmax;
	ras->y0		=
	ras->ymin	= ymin;
	ras->y1		=
	ras->ymax	= ymax;
	ras->xstep	= (xmax-xmin)/(ncol - 2 * mcol - 1);
	ras->ystep	= (ymax-ymin)/(nrow - 2 * mrow - 1);
	/* Rprintf("xstep,ystep = %lf,%lf\n", ras->xstep,ras->ystep);  */
}
