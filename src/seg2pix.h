/*
  seg2pix.h

  Code template for seg2pix.c

  $Revision: 1.2 $ $Date: 2015/01/08 10:57:20 $

  Macros:
  FNAME   name of function
  SUMUP   #defined if crossings should be counted (weights summed)
 
  V       matrix index macro (in seg2pix.c)
  DEBUG   debug if #defined

*/

#undef INCREMENT
#undef ZERO

#ifdef SUMUP
#define ZERO (double) 0.0
#define INCREMENT(I,J) V(I,J) += wi
#else
#define ZERO 0
#define INCREMENT(I,J) V(I,J) = 1
#endif


void FNAME(ns,x0,y0,x1,y1,
#ifdef SUMUP
	   w,
#endif
	   nx,ny,out)
     int *ns;  /* number of segments */
     double *x0,*y0,*x1,*y1; /* coordinates of segment endpoints */
     int *nx, *ny;  /* dimensions of pixel array (columns, rows) */
#ifdef SUMUP
     double *w; /* weights attached to segments */
     double *out; /* output totals */
#else 
     int *out;     /* output indicators */
#endif
{
  int Ns, Nx, Ny, i, j, k, m, m0, m1, mmin, mmax, maxchunk;
  double x0i, x1i, y0i, y1i, dx, dy;
  double leni;
  double xleft, yleft, xright, yright, slope;
  double xstart, ystart, xfinish, yfinish;
  int mleft, mright, kstart, kfinish, kmin, kmax;
#ifdef SUMUP
  double wi;
#endif

  Ns = *ns;
  Nx = *nx;
  Ny = *ny;
  
  for(k = 0; k < Ny - 1; k++) 
    for(j = 0; j < Nx - 1; j++) 
      V(k, j) = ZERO;

  OUTERCHUNKLOOP(i, Ns, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Ns, maxchunk, 8196) {
      x0i = x0[i];
      y0i = y0[i];
      x1i = x1[i];
      y1i = y1[i];   
#ifdef SUMUP
      wi = w[i];
#endif
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
	INCREMENT(j, k);
      } else if(floor(x1i) == floor(x0i) && floor(y1i) == floor(y0i)) { 
	/* contained in one cell */
#ifdef DEBUG
	Rprintf("contained in one cell\n");
#endif
	k = clamp((int) floor(x0i), 0, Nx-1);
	j = clamp((int) floor(y0i), 0, Ny-1);
	INCREMENT(j, k);
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
	  INCREMENT(j,k);
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
	  INCREMENT(j,k);
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
	    INCREMENT(k, m);
	}
      } /* end of if-else */
    }
  }
#ifdef DEBUG
    Rprintf("done\n");
#endif
}

