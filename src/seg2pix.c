#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "chunkloop.h"

#undef DEBUG

/*

  seg2pix.c

  Discretise line segment on pixel grid

  seg2pixI      pixel value is indicator = 1 if any line crosses pixel

  seg2pixL      pixel value is total (weighted) length of lines inside pixel

  (rescale R data so that pixels are integer)
  pixels numbered 0, ..., nx-1 and 0, ..., ny-1
  with boundaries at x=0, x=nx, y=0, y=ny.

*/

#define V(I,J) out[(I) + (J) * (Ny)]

int clamp(k, n0, n1) 
     int k, n0, n1;
{
  int m;
  m = k;
  if(m < n0) m = n0; 
  if(m > n1) m = n1;
  return(m);
}

void seg2pixI(ns,x0,y0,x1,y1,nx,ny,out)
     int *ns;  /* number of segments */
     double *x0,*y0,*x1,*y1; /* coordinates of segment endpoints */
     int *nx, *ny;  /* dimensions of pixel array (columns, rows) */
     int *out;     
{
  int Ns, Nx, Ny, i, j, k, m, m0, m1, mmin, mmax, maxchunk;
  double x0i, x1i, y0i, y1i, dx, dy;
  double leni;
  double xleft, yleft, xright, yright, slope;
  double xstart, ystart, xfinish, yfinish;
  int mleft, mright, kstart, kfinish, kmin, kmax;

  Ns = *ns;
  Nx = *nx;
  Ny = *ny;
  
  for(k = 0; k < Ny - 1; k++)
    for(j = 0; j < Nx - 1; j++)
      V(k, j) = 0;

  OUTERCHUNKLOOP(i, Ns, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Ns, maxchunk, 8196) {
      x0i = x0[i];
      y0i = y0[i];
      x1i = x1[i];
      y1i = y1[i];   
      dx = x1i - x0i;
      dy = y1i - y0i;
      leni = hypot(dx, dy);
#ifdef DEBUG
      Rprintf("(%lf, %lf) to (%lf, %lf)\n",
	      x0i, y0i, x1i, y1i);
#endif
      if(leni < 0.001) { /* tiny segment */
#ifdef DEBUG
	Rprintf("tiny\n");
#endif
	k = clamp((int) floor(x0i), 0, Nx-1);
	j = clamp((int) floor(y0i), 0, Ny-1);
	V(j,k) = 1;
      } else if(floor(x1i) == floor(x0i) && floor(y1i) == floor(y0i)) { 
	/* contained in one cell */
#ifdef DEBUG
	Rprintf("contained in one cell\n");
#endif
	k = clamp((int) floor(x0i), 0, Nx-1);
	j = clamp((int) floor(y0i), 0, Ny-1);
	V(j,k) = 1;
      } else if(floor(y1i) == floor(y0i)) { /* horizontal */
#ifdef DEBUG
	Rprintf("horizontal\n");
#endif
	j = clamp((int) floor(y1i), 0, Ny-1);
	m0 = clamp((int) floor(x0i), 0, Nx-1);
	m1 = clamp((int) floor(x1i), 0, Nx-1);
	mmin = (m0 < m1) ? m0: m1;
	mmax = (m0 < m1) ? m1: m0;
#ifdef DEBUG
	Rprintf("row %d: columns [%d, %d]\n", j, mmin, mmax);
#endif
	for(k = mmin; k <= mmax; k++) 
	  V(j,k) = 1;
      } else if(floor(x1i) == floor(x0i)) { /* vertical */
#ifdef DEBUG
	Rprintf("vertical\n");
#endif
	k = clamp((int) floor(x1i), 0, Nx-1);
	m0 = clamp((int) floor(y0i), 0, Ny-1);
	m1 = clamp((int) floor(y1i), 0, Ny-1);
	mmin = (m0 < m1) ? m0: m1;
	mmax = (m0 < m1) ? m1: m0;
#ifdef DEBUG
	Rprintf("column %d: rows [%d, %d]\n", k, mmin, mmax);
#endif
	for(j = mmin; j <= mmax; j++) 
	  V(j,k) = 1;
      } else { /* general case */
#ifdef DEBUG
	Rprintf("general\n");
#endif
	if(x1i > x0i) {
	  xleft = x0i;
	  yleft = y0i;
	  xright = x1i;
	  yright = y1i;
	} else {
	  xleft = x1i;
	  yleft = y1i;
	  xright = x0i;
	  yright = y0i;
	}
	slope = (yright - yleft)/(xright - xleft);
	mleft = clamp((int) floor(xleft), 0, Nx-1);
	mright = clamp((int) floor(xright), 0, Nx-1); 
#ifdef DEBUG
	Rprintf("column range [%d, %d]\n", mleft, mright);
#endif
	/* treat each vertical slice */
	for(m = mleft; m <= mright; m++) {
	  if(m == mleft) {
	    xstart = xleft;
	    ystart = yleft;
	  } else {
	    xstart = m;
	    ystart = yleft + slope * (xstart - xleft);
	  }
	  if(m == mright) {
	    xfinish = xright;
	    yfinish = yright;
	  } else {
	    xfinish = m+1;
	    yfinish = yleft + slope * (xfinish - xleft);
	  }
	  kstart = clamp((int) floor(ystart), 0, Ny-1);
	  kfinish = clamp((int) floor(yfinish), 0, Ny-1);
	  kmin = (kstart < kfinish) ? kstart : kfinish;
	  kmax = (kstart < kfinish) ? kfinish : kstart;
#ifdef DEBUG
	  Rprintf("column %d: rows [%d, %d]\n", m, kmin, kmax);
#endif
	  for(k = kmin; k <= kmax; k++)
	    V(k, m) = 1;
	}
      } /* end of if-else */
    }
  }
#ifdef DEBUG
    Rprintf("done\n");
#endif
}


void seg2pixL(ns,x0,y0,x1,y1,weights,pixwidth,pixheight,nx,ny,out)
     int *ns;
     double *x0,*y0,*x1,*y1,*weights; /* segment coordinates and weights */
     double *pixwidth, *pixheight;  /* original pixel dimensions */
     int *nx, *ny;
     double *out;  /* output matrix */
{
  int Ns, Nx, Ny, i, j, k, m, mmin, mmax, maxchunk;
  double x0i, x1i, y0i, y1i;
  double leni;
  double xleft, yleft, xright, yright, slope, scalesecant;
  double xlow, xhigh, ylow, yhigh, invslope, scalecosecant;
  double xstart, ystart, xfinish, yfinish; 
  double xxx0, xxx1, yyy0, yyy1;
  int mleft, mright, kstart, kfinish, kmin, kmax;
  double pwidth, pheight, pwidth2, pheight2;
  double wti; 

  Ns = *ns;
  Nx = *nx;
  Ny = *ny;

  /* 
     one scaled x unit = 'pwidth' original x units
     one scaled y unit = 'pheight' original y units
  */
	 
  pwidth = *pixwidth;
  pheight = *pixheight;
  pwidth2 = pwidth * pwidth;
  pheight2 = pheight * pheight;

  /* zero the matrix */

  for(k = 0; k < Ny - 1; k++)
    for(j = 0; j < Nx - 1; j++)
      V(k, j) = 0;

  OUTERCHUNKLOOP(i, Ns, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Ns, maxchunk, 8196) {
      x0i = x0[i];
      y0i = y0[i];
      x1i = x1[i];
      y1i = y1[i];   
      wti = weights[i];
      leni = sqrt(pwidth2 * pow(x1i - x0i, 2) + pheight2 * pow(y1i-y0i, 2));
#ifdef DEBUG
      Rprintf("(%lf, %lf) to (%lf, %lf), length %lf\n",
	      x0i, y0i, x1i, y1i, leni);
#endif
      if(leni < 0.001) { /* tiny segment */
#ifdef DEBUG
	Rprintf("tiny\n");
#endif
	k = clamp((int) floor(x0i), 0, Nx-1);
	j = clamp((int) floor(y0i), 0, Ny-1);
	V(j,k) += wti * leni;
      } else if(floor(x1i) == floor(x0i) && floor(y1i) == floor(y0i)) { 
	/* contained in one cell */
#ifdef DEBUG
	Rprintf("contained in one cell\n");
#endif
	k = clamp((int) floor(x0i), 0, Nx-1);
	j = clamp((int) floor(y0i), 0, Ny-1);
	V(j,k) += wti * leni;
      } else if(floor(y1i) == floor(y0i)) { /* horizontal */
#ifdef DEBUG
	Rprintf("horizontal\n");
#endif
	j = clamp((int) floor(y1i), 0, Ny-1);
	if(x1i > x0i) {
	  xleft = x0i;
	  yleft = y0i;
	  xright = x1i;
	  yright = y1i;
	} else {
	  xleft = x1i;
	  yleft = y1i;
	  xright = x0i;
	  yright = y0i;
	}
	mmin = clamp((int) floor(xleft), 0, Nx-1);
	mmax = clamp((int) floor(xright), 0, Nx-1);
	slope = (yright - yleft)/(xright - xleft);
	scalesecant = wti * sqrt(pwidth2 + slope * slope * pheight2);
	/* 
	   For this slope, one scaled x unit means
	   'pwidth' original x units and
	   slope * pheight original y units
	   i.e. line length sqrt(pwidth^2 + slope^2 * pheight^2)
	 
	*/
	for(k = mmin; k <= mmax; k++) {
	  xstart = (k == mmin) ? xleft : k;
	  xfinish = (k == mmax) ? xright : (k+1);
	  V(j,k) += (xfinish - xstart) * scalesecant;
	}
      } else if(floor(x1i) == floor(x0i)) { /* vertical */
#ifdef DEBUG
	Rprintf("vertical\n");
#endif
	k = clamp((int) floor(x1i), 0, Nx-1);
	if(y1i > y0i) {
	  xlow = x0i;
	  ylow = y0i;
	  xhigh = x1i;
	  yhigh = y1i;
	} else {
	  xlow = x1i;
	  ylow = y1i;
	  xhigh = x0i;
	  yhigh = y0i;
	}
	mmin = clamp((int) floor(ylow), 0, Ny-1);
	mmax = clamp((int) floor(yhigh), 0, Ny-1);
	invslope = (xhigh - xlow)/(yhigh - ylow);
	scalecosecant = wti * sqrt(pheight2 + invslope * invslope * pwidth2);
#ifdef DEBUG
	Rprintf("i = %d\n", i);
	Rprintf("inverse slope = %lf\n", invslope);
	Rprintf("scaled cosecant = %lf\n", scalecosecant);
#endif
	/* 
	   For this slope, one scaled y unit means
	   'pheight' original y units and
	   invslope * pwidth original x units
	   i.e. line length sqrt(pheight^2 + invslope^2 * pwidth^2)
	 
	*/
	for(j = mmin; j <= mmax; j++) {
	  ystart = (j == mmin)? ylow : j;
	  yfinish = (j == mmax)? yhigh : (j+1);
	  V(j,k) += (yfinish - ystart) * scalecosecant;
	}
      } else { /* general case */
#ifdef DEBUG
	Rprintf("general\n");
#endif
	if(x1i > x0i) {
	  xleft = x0i;
	  yleft = y0i;
	  xright = x1i;
	  yright = y1i;
	} else {
	  xleft = x1i;
	  yleft = y1i;
	  xright = x0i;
	  yright = y0i;
	}
	slope = (yright - yleft)/(xright - xleft);
	mleft = clamp((int) floor(xleft), 0, Nx-1);
	mright = clamp((int) floor(xright), 0, Nx-1); 
#ifdef DEBUG
	Rprintf("column range [%d, %d]\n", mleft, mright);
#endif
	/* treat each vertical slice */
	for(m = mleft; m <= mright; m++) {
	  if(m == mleft) {
	    xstart = xleft;
	    ystart = yleft;
	  } else {
	    xstart = m;
	    ystart = yleft + slope * (xstart - xleft);
	  }
	  if(m == mright) {
	    xfinish = xright;
	    yfinish = yright;
	  } else {
	    xfinish = m+1;
	    yfinish = yleft + slope * (xfinish - xleft);
	  }
	  kstart = clamp((int) floor(ystart), 0, Ny-1);
	  kfinish = clamp((int) floor(yfinish), 0, Ny-1);
	  if(ystart < yfinish) {
	    kmin = kstart;
	    kmax = kfinish;
	    ylow = ystart;
	    yhigh = yfinish;
	  } else {
	    kmin = kfinish;
	    kmax = kstart;
	    ylow = yfinish;
	    yhigh = ystart;
	  }
#ifdef DEBUG
	  Rprintf("column %d: rows [%d, %d]\n", m, kmin, kmax);
#endif
	  for(k = kmin; k <= kmax; k++) { 
	    yyy0 = (k == kmin) ? ylow : k;
	    yyy1 = (k == kmax) ? yhigh : (k+1);
	    xxx0 = xstart + (yyy0 - ystart)/slope;
	    xxx1 = xstart + (yyy1 - ystart)/slope;
	    V(k, m) += wti * sqrt(pow(yyy1 - yyy0, 2) * pheight2 + 
				  pow(xxx1 - xxx0, 2) * pwidth2);
	  }
	}
      }
    }
  }
#ifdef DEBUG
  Rprintf("done.\n");
#endif
}



