/*
       distseg.c

       Distances from point pattern to line segment pattern
       Distance transform of a line segment pattern

       nndist2segs: minimum distance from point to any line segment
       prdist2segs: pairwise distances from each point to each line segment

       $Revision: 1.9 $ $Date: 2012/03/27 05:38:51 $

       Author: Adrian Baddeley

*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "chunkloop.h"

void
nndist2segs(xp, yp, npoints, x0, y0, x1, y1, nsegments, epsilon, dist2, index)
     /* input */
     double	*xp, *yp;		/* point/pixel coordinates */
     int	*npoints;
     double	*x0, *y0, *x1, *y1;	/* line segment endpoints */
     int	*nsegments;
     double     *epsilon;               /* tolerance for short segments */
     /* output */
     double	*dist2;		        /* squared distance from pixel 
                                        to nearest line segment 
					  INITIALISED TO LARGE VALUE 
					*/
     int	*index;		        /* which line segment is closest */
{
  int	i,j, np, nseg, maxchunk;
  double dx,dy,leng,co,si;  /* parameters of segment */
  double xdif0,ydif0,xdif1,ydif1,xpr,ypr; /* vectors */
  double dsq0,dsq1,dsq,dsqperp; /* squared distances */
  double eps;

  np   = *npoints;
  nseg = *nsegments;
  eps  = *epsilon;

  OUTERCHUNKLOOP(j, nseg, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nseg, maxchunk, 16384) {
      dx = x1[j] - x0[j];
      dy = y1[j] - y0[j];
      leng = hypot(dx, dy);
      if(leng > eps) {
	/* normal case */
	co = dx/leng;
	si = dy/leng;
	for(i = 0; i < np; i++) {
	  /* vectors from pixel to segment endpoints */
	  xdif0 =  xp[i] - x0[j];
	  ydif0 =  yp[i] - y0[j];
	  xdif1 =  xp[i] - x1[j];
	  ydif1 =  yp[i] - y1[j];
	  /* squared distances to segment endpoints */
	  dsq0 = xdif0 * xdif0 + ydif0 * ydif0;
	  dsq1 = xdif1 * xdif1 + ydif1 * ydif1;
	  dsq = (dsq0 < dsq1) ? dsq0 : dsq1;
	  /* rotate pixel around 1st endpoint of segment
	     so that line segment lies in x axis */
	  xpr = xdif0 * co + ydif0 * si;
	  ypr = -xdif0 * si + ydif0 * co;
	  /* perpendicular distance applies only in perpendicular region */
	  if(xpr >= 0.0 && xpr <= leng) {
	    dsqperp = ypr * ypr;
	    if(dsqperp < dsq) dsq = dsqperp;
	  }
	  if(dist2[i] > dsq) {
	    dist2[i] = dsq;
	    index[i] = j;
	  }
	}
      } else {
	/* short segment - use endpoints only */
	for(i = 0; i < np; i++) {
	  /* vectors from pixel to segment endpoints */
	  xdif0 =  xp[i] - x0[j];
	  ydif0 =  yp[i] - y0[j];
	  xdif1 =  xp[i] - x1[j];
	  ydif1 =  yp[i] - y1[j];
	  /* squared distances to segment endpoints */
	  dsq0 = xdif0 * xdif0 + ydif0 * ydif0;
	  dsq1 = xdif1 * xdif1 + ydif1 * ydif1;
	  dsq = (dsq0 < dsq1) ? dsq0 : dsq1;
	  if(dist2[i] > dsq) {
	    dist2[i] = dsq;
	    index[i] = j;
	  }
	}
      }
    }
  }
}

void
prdist2segs(xp, yp, npoints, x0, y0, x1, y1, nsegments, epsilon, dist2)
     /* input */
     double	*xp, *yp;		/* point/pixel coordinates */
     int	*npoints;
     double	*x0, *y0, *x1, *y1;	/* line segment endpoints */
     int	*nsegments;
     double     *epsilon;               /* tolerance for short segments */
     /* output */
     double	*dist2;		        /* squared distances from each pixel 
                                        to each line segment */
{
  int	i,j, np, nseg, maxchunk;
  double dx,dy,leng,co,si;  /* parameters of segment */
  double xdif0,ydif0,xdif1,ydif1,xpr,ypr; /* vectors */
  double dsq0,dsq1,dsq,dsqperp; /* squared distances */
  double eps;

  np   = *npoints;
  nseg = *nsegments;
  eps  = *epsilon;

  OUTERCHUNKLOOP(j, nseg, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nseg, maxchunk, 16384) {
      dx = x1[j] - x0[j];
      dy = y1[j] - y0[j];
      leng = hypot(dx, dy);
      if(leng > eps) {
	/* normal case */
	co = dx/leng;
	si = dy/leng;
	for(i = 0; i < np; i++) {
	  /* vectors from pixel to segment endpoints */
	  xdif0 =  xp[i] - x0[j];
	  ydif0 =  yp[i] - y0[j];
	  xdif1 =  xp[i] - x1[j];
	  ydif1 =  yp[i] - y1[j];
	  /* squared distances to segment endpoints */
	  dsq0 = xdif0 * xdif0 + ydif0 * ydif0;
	  dsq1 = xdif1 * xdif1 + ydif1 * ydif1;
	  dsq = (dsq0 < dsq1) ? dsq0 : dsq1;
	  /* rotate pixel around 1st endpoint of segment
	     so that line segment lies in x axis */
	  xpr = xdif0 * co + ydif0 * si;
	  ypr = -xdif0 * si + ydif0 * co;
	  /* perpendicular distance applies only in perpendicular region */
	  if(xpr >= 0.0 && xpr <= leng) {
	    dsqperp = ypr * ypr;
	    if(dsqperp < dsq) dsq = dsqperp;
	  }
	  dist2[i + j * np] = dsq;
	}
      } else {
	/* short segment */
	for(i = 0; i < np; i++) {
	  /* vectors from pixel to segment endpoints */
	  xdif0 =  xp[i] - x0[j];
	  ydif0 =  yp[i] - y0[j];
	  xdif1 =  xp[i] - x1[j];
	  ydif1 =  yp[i] - y1[j];
	  /* squared distances to segment endpoints */
	  dsq0 = xdif0 * xdif0 + ydif0 * ydif0;
	  dsq1 = xdif1 * xdif1 + ydif1 * ydif1;
	  dsq = (dsq0 < dsq1) ? dsq0 : dsq1;
	  dist2[i + j * np] = dsq;
	}
      }
    }
  }	
}

