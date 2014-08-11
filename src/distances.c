/*

  distances.c

  Distances between pairs of points

  $Revision: 1.27 $     $Date: 2012/03/27 04:22:13 $

  Cpairdist      Pairwise distances
  Cpair2dist     Pairwise distances squared
  CpairPdist     Pairwise distances with periodic correction
  CpairP2dist    Pairwise distances squared, with periodic correction

  Ccrossdist     Pairwise distances for two sets of points
  Ccross2dist    Pairwise distances squared, for two sets of points
  CcrossPdist    Pairwise distances for two sets of points, periodic correction

  Cmatchxy       Find matches between two sets of points   

 */

#include <math.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

double sqrt();

void Cpairdist(n, x, y, d)
     /* inputs */
     int *n;
     double *x, *y;
     /* output */
     double *d;
{ 
  int i, j, npoints, maxchunk; 
  double *dp;
  double xi, yi, dx, dy, dist;

  npoints = *n;

  /* set d[0,0] = 0 */
  *d = 0.0;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {
      xi = x[i];
      yi = y[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dist = sqrt( dx * dx + dy * dy ); 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
  }
}

/* squared distances */

void Cpair2dist(n, x, y, d)
     /* inputs */
     int *n;
     double *x, *y;
     /* output */
     double *d;
{ 
  int i, j, npoints, maxchunk; 
  double *dp;
  double xi, yi, dx, dy, dist;

  npoints = *n;

  /* set d[0,0] = 0 */
  *d = 0.0;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {
      xi = x[i];
      yi = y[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dist = dx * dx + dy * dy; 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
  }
}

void Ccrossdist(nfrom, xfrom, yfrom, nto, xto, yto, d)
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *xto, *yto;
     /* output */
     double *d;
{ 
  int i, j, nf, nt, maxchunk; 
  double *dptr;
  double xj, yj, dx, dy;

  nf = *nfrom;
  nt = *nto;

  dptr = d;

  OUTERCHUNKLOOP(j, nt, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nt, maxchunk, 16384) {
      xj = xto[j];
      yj = yto[j];
      for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	*dptr = sqrt( dx * dx + dy * dy ); 
      }
    }
  }
}


/* squared distances */

void Ccross2dist(nfrom, xfrom, yfrom, nto, xto, yto, d)
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *xto, *yto;
     /* output */
     double *d;
{ 
  int i, j, nf, nt, maxchunk; 
  double *dptr;
  double xj, yj, dx, dy;

  nf = *nfrom;
  nt = *nto;

  dptr = d;

  OUTERCHUNKLOOP(j, nt, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nt, maxchunk, 16384) {
      xj = xto[j];
      yj = yto[j];
      for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	*dptr = dx * dx + dy * dy; 
      }
    }
  }
}


/* distances with periodic correction */

void CpairPdist(n, x, y, xwidth, yheight, d)
     /* inputs */
     int *n;
     double *x, *y, *xwidth, *yheight;
     /* output */
     double *d;
{ 
  int i, j, npoints, maxchunk; 
  double *dp;
  double xi, yi, dx, dy, dx2, dy2, dx2p, dy2p, dist, wide, high;

  npoints = *n;
  wide = *xwidth;
  high = *yheight;

  /* set d[0,0] = 0 */
  *d = 0.0;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {
      xi = x[i];
      yi = y[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dist = sqrt( dx2p + dy2p ); 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
  }
}

/* same function without the sqrt */

void CpairP2dist(n, x, y, xwidth, yheight, d)
     /* inputs */
     int *n;
     double *x, *y, *xwidth, *yheight;
     /* output */
     double *d;
{ 
  int i, j, npoints, maxchunk; 
  double *dp;
  double xi, yi, dx, dy, dx2, dy2, dx2p, dy2p, dist, wide, high;

  npoints = *n;
  wide = *xwidth;
  high = *yheight;

  /* set d[0,0] = 0 */
  *d = 0.0;

  OUTERCHUNKLOOP(i, npoints, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, npoints, maxchunk, 16384) {
      xi = x[i];
      yi = y[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dist = dx2p + dy2p; 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
  }
}

void CcrossPdist(nfrom, xfrom, yfrom, nto, xto, yto, xwidth, yheight, d)
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *xto, *yto, *xwidth, *yheight;
     /* output */
     double *d;
{ 
  int i, j, nf, nt, maxchunk; 
  double *dptr;
  double xj, yj, dx, dy, dx2, dy2, dx2p, dy2p, wide, high;

  nf = *nfrom;
  nt = *nto;
  wide = *xwidth;
  high = *yheight;

  dptr = d;

  OUTERCHUNKLOOP(j, nt, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nt, maxchunk, 16384) {
      xj = xto[j];
      yj = yto[j];
      for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	dx2p = dx * dx;
	dy2p = dy * dy;
	dx2 = (dx - wide) * (dx - wide);
	dy2 = (dy - high) * (dy - high);
	if(dx2 < dx2p) dx2p = dx2;
	if(dy2 < dy2p) dy2p = dy2;
	dx2 = (dx + wide) * (dx + wide);
	dy2 = (dy + high) * (dy + high);
	if(dx2 < dx2p) dx2p = dx2;
	if(dy2 < dy2p) dy2p = dy2;
	*dptr = sqrt( dx2p + dy2p ); 
      }
    }
  }
}

/*

  matchxy

  Find matches between two lists of points

 */

void Cmatchxy(na, xa, ya, nb, xb, yb, match)
     /* inputs */
     int *na, *nb;
     double *xa, *ya, *xb, *yb;
     /* output */
     int *match; 
     /* match[i] = j+1 if xb[j], yb[j] matches xa[i], ya[i] */
     /* match[i] = 0   if no such point matches xa[i], ya[i] */
{ 
  int i, j, Na, Nb, maxchunk; 
  double xai, yai;

  Na = *na;
  Nb = *nb;

  OUTERCHUNKLOOP(i, Na, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Na, maxchunk, 16384) {
      xai = xa[i];
      yai = ya[i];
      match[i] = 0;
      for (j=0; j < Nb; j++) {
	if(xai == xb[j] && yai == yb[j]) {
	  match[i] = j+1;
	  break;
	}
      }
    }
  }
}
