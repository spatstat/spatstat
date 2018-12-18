/*
       metricPdist.h

       Distance transform of a discrete binary image
       using a general metric
       
       Code template which is #included several times in metricPdist.c

       $Revision: 1.3 $ $Date: 2018/12/18 02:43:11 $

       Uses the following definitions
       FNAME          Function name (called from R)
       MARGLIST       List of function arguments specifying the metric
       MARGDECLARE    Declarations of function arguments specifying the metric
       MTEMPDECLARE   Declaration and initialisation of variables for metric
       METRIC         Expression for calculating the metric (x1,y1,x2,y2)

       Also uses definitions from raster.h and metricPdist.c

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

void FNAME(xmin, ymin, xmax, ymax,
	   nr, nc, mr, mc, 
	   inp,
	   MARGLIST,
	   npasses,
	   distances, rows, cols
	   )
	double *xmin, *ymin, *xmax, *ymax;  	  /* x, y dimensions */
	int *nr, *nc;	 	                  /* raster dimensions
				                     EXCLUDING margins */
	int *mr, *mc;                             /* margins */
	int   *inp;            /* input:  binary image */
        MARGDECLARE;
	int *npasses;          /* number of passes over raster */
	double *distances;     /* output: distance to nearest point */
	int   *rows;	       /* output: row of nearest point (start= 0) */
	int   *cols;	       /* output: column of nearest point (start = 0) */
	/* all images must have identical dimensions 
	   including a margin of 1 on each side */
{
	Raster data, dist, row, col;
	int mrow, mcol, nrow, ncol;

	int	j,k;
	double	d, x, y;
	int	r, c;
	int     Npass, ipass;
	double	dnew, diam, dd, huge;
	double  Xmin, Ymin, Xmax, Ymax;
	
	/* declare any variables used for the metric */
	MTEMPDECLARE;

        Xmin = *xmin;
        Xmax = *xmax;
        Ymin = *ymin;
        Ymax = *ymax;
	mrow = *mr;
	mcol = *mc;

	Npass = *npasses;

	/* Determine diameter of window. 
	   (must be achieved as distance between two of the vertices)
	*/
	/* diagonals */
	diam = METRIC(Xmin,Ymin,Xmax,Ymax); 
	dd   = METRIC(Xmin,Ymax,Xmax,Ymin);
	if(dd > diam) diam = dd;
	dd = METRIC(Xmin,Ymin,Xmin,Ymax);
	if(dd > diam) diam = dd;
	/* horizontals */
	dd = METRIC(Xmin,Ymin,Xmax,Ymin);
	if(dd > diam) diam = dd;
	dd = METRIC(Xmin,Ymax,Xmax,Ymax);
	if(dd > diam) diam = dd;
	/* verticals */
	dd = METRIC(Xmin,Ymin,Xmin,Ymax);
	if(dd > diam) diam = dd;
	dd = METRIC(Xmax,Ymin,Xmax,Ymax);
	if(dd > diam) diam = dd;

	/* create raster structures */

	/* full dimensions */
	nrow = *nr + 2 * mrow;
	ncol = *nc + 2 * mcol;

	shape_raster( &data, (void *) inp,
		      Xmin, Ymin, Xmax, Ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &dist, (void *) distances,
		      Xmin, Ymin, Xmax, Ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &row, (void *) rows,
		      Xmin, Ymin, Xmax, Ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &col, (void *) cols,
		      Xmin, Ymin, Xmax, Ymax,
		      nrow, ncol, mrow, mcol);

	/* initialise arrays */
	
	Clear(row,int,UNDEFINED)
	Clear(col,int,UNDEFINED)
	huge = 2.0 * diam;
	Clear(dist,double,huge)

	/* 
	   if input pixel is TRUE, 
	   set distance to 0 and make pixel point to itself 
	*/
	for(j = data.rmin; j <= data.rmax; j++)
	for(k = data.cmin; k <= data.cmax; k++) 
	  if(Entry(data, j, k, int) != 0) {
	      Entry(dist, j, k, double) = 0.0;
	      Entry(row,  j, k, int)   = j;
	      Entry(col,  j, k, int)   = k;
	  }

	/* how to update the distance values */

#undef GETVALUES	
#define GETVALUES(ROW,COL) \
	x = Xpos(data, COL); \
	y = Ypos(data, ROW); \
	d = Entry(dist,ROW,COL,double); 

#undef COMPARE
#define COMPARE(ROW,COL,RR,CC) \
	r = Entry(row,RR,CC,int); \
	c = Entry(col,RR,CC,int); \
	if(Is_Defined(r) && Is_Defined(c) \
	   && Entry(dist,RR,CC,double) < d) { \
	     dnew = METRIC(x, y, Xpos(data,c), Ypos(data,r)); \
	     if(dnew < d) { \
		Entry(row,ROW,COL,int) = r; \
		Entry(col,ROW,COL,int) = c; \
		Entry(dist,ROW,COL,double) = dnew; \
		d = dnew; \
	     } \
	}

	for(ipass = 0; ipass < Npass; ipass++) {
	  
	/* forward pass */

	for(j = data.rmin; j <= data.rmax; j++)
	for(k = data.cmin; k <= data.cmax; k++) {
	        GETVALUES(j, k)
		COMPARE(j,k, j-1,k-1)
		COMPARE(j,k, j-1,  k)
		COMPARE(j,k, j-1,k+1)
		COMPARE(j,k, j,  k-1)
		  }

	/* backward pass */

	for(j = data.rmax; j >= data.rmin; j--) 
	for(k = data.cmax; k >= data.cmin; k--) {
	        GETVALUES(j, k)
		COMPARE(j,k, j+1,k+1)
		COMPARE(j,k, j+1,  k)
		COMPARE(j,k, j+1,k-1)
		COMPARE(j,k, j,  k+1)
		  }

	}
}

