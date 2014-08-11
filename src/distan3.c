/*

  distan3.c

  Distances between pairs of 3D points

  $Revision: 1.1 $     $Date: 2010/01/05 05:03:31 $

  D3pairdist      Pairwise distances
  D3pair2dist     Pairwise distances squared
  D3pairPdist     Pairwise distances with periodic correction
  D3pairP2dist    Pairwise distances squared, with periodic correction

  D3crossdist     Pairwise distances for two sets of points
  D3cross2dist    Pairwise distances squared, for two sets of points
  D3crossPdist    Pairwise distances for two sets of points, periodic correction

  matchxyz       Find matches between two sets of points   

 */

#include <math.h>
/* #include <stdio.h> */

double sqrt();

void D3pairdist(n, x, y, z, d)
     /* inputs */
     int *n;
     double *x, *y, *z;
     /* output */
     double *d;
{ 
  int i, j, npoints; 
  double *dp;
  double xi, yi, zi, dx, dy, dz, dist;

  npoints = *n;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dz = z[j] - zi;
	  dist = sqrt( dx * dx + dy * dy + dz * dz ); 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

/* squared distances */

void D3pair2dist(n, x, y, z, d)
     /* inputs */
     int *n;
     double *x, *y, *z;
     /* output */
     double *d;
{ 
  int i, j, npoints; 
  double *dp;
  double xi, yi, zi, dx, dy, dz, dist;

  npoints = *n;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dz = z[j] - zi;
	  dist = dx * dx + dy * dy + dz * dz; 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

void D3crossdist(nfrom, xfrom, yfrom, zfrom, nto, xto, yto, zto, d)
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *zfrom, *xto, *yto, *zto;
     /* output */
     double *d;
{ 
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, zj, dx, dy, dz;

  nf = *nfrom;
  nt = *nto;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    zj = zto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	dz = zj - zfrom[i];
	*dptr = sqrt( dx * dx + dy * dy + dz * dz ); 
    }
  }
}

/* squared distances */

void D3cross2dist(nfrom, xfrom, yfrom, zfrom, nto, xto, yto, zto, d)
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *zfrom, *xto, *yto, *zto;
     /* output */
     double *d;
{ 
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, zj, dx, dy, dz;

  nf = *nfrom;
  nt = *nto;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    zj = zto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	dz = zj - zfrom[i];
	*dptr = dx * dx + dy * dy + dz * dz; 
    }
  }
}



/* distances with periodic correction */

void D3pairPdist(n, x, y, z, xwidth, yheight, zdepth, d)
     /* inputs */
     int *n;
     double *x, *y, *z, *xwidth, *yheight, *zdepth;
     /* output */
     double *d;
{ 
  int i, j, npoints; 
  double *dp;
  double xi, yi, zi, dx, dy, dz, dx2, dy2, dz2, dx2p, dy2p, dz2p, dist, wide, high, deep;

  npoints = *n;
  wide = *xwidth;
  high = *yheight;
  deep = *zdepth;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dz = z[j] - zi;
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dz2p = dz * dz;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  dz2 = (dz - deep) * (dz - deep);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  if(dz2 < dz2p) dz2p = dz2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  dz2 = (dz + deep) * (dz + deep);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  if(dz2 < dz2p) dz2p = dz2;
	  dist = sqrt( dx2p + dy2p + dz2p ); 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

/* same function without the sqrt */

void D3pairP2dist(n, x, y, z, xwidth, yheight, zdepth, d)
     /* inputs */
     int *n;
     double *x, *y, *z, *xwidth, *yheight, *zdepth;
     /* output */
     double *d;
{ 
  int i, j, npoints; 
  double *dp;
  double xi, yi, zi, dx, dy, dz, dx2, dy2, dz2, dx2p, dy2p, dz2p, dist, wide, high, deep;

  npoints = *n;
  wide = *xwidth;
  high = *yheight;
  deep = *zdepth;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dz = z[j] - zi;
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dz2p = dz * dz;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  dz2 = (dz - deep) * (dz - deep);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  if(dz2 < dz2p) dz2p = dz2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  dz2 = (dz + deep) * (dz + deep);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  if(dz2 < dz2p) dz2p = dz2;
	  dist = dx2p + dy2p + dz2p; 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

void D3crossPdist(nfrom, xfrom, yfrom, zfrom, nto, xto, yto, zto, xwidth, yheight, zdepth, d)
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *zfrom, *xto, *yto, *zto, *xwidth, *yheight, *zdepth;
     /* output */
     double *d;
{ 
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, zj, dx, dy, dz, dx2, dy2, dz2, dx2p, dy2p, dz2p, wide, high, deep;

  nf = *nfrom;
  nt = *nto;
  wide = *xwidth;
  high = *yheight;
  deep = *zdepth;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    zj = zto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	dz = zj - zfrom[i];
	dx2p = dx * dx;
	dy2p = dy * dy;
	dz2p = dz * dz;
	dx2 = (dx - wide) * (dx - wide);
	dy2 = (dy - high) * (dy - high);
	dz2 = (dz - deep) * (dz - deep);
	if(dx2 < dx2p) dx2p = dx2;
	if(dy2 < dy2p) dy2p = dy2;
	if(dz2 < dz2p) dz2p = dz2;
	dx2 = (dx + wide) * (dx + wide);
	dy2 = (dy + high) * (dy + high);
	dz2 = (dy + deep) * (dz + deep);
	if(dx2 < dx2p) dx2p = dx2;
	if(dy2 < dy2p) dy2p = dy2;
	if(dz2 < dz2p) dz2p = dz2;
	*dptr = sqrt( dx2p + dy2p + dz2p ); 
    }
  }
}

/*

  matchxyz

  Find matches between two lists of points

 */

void matchxyz(na, xa, ya, za, nb, xb, yb, zb, match)
     /* inputs */
     int *na, *nb;
     double *xa, *ya, *za, *xb, *yb, *zb;
     /* output */
     int *match; 
{ 
  int i, j, Na, Nb; 
  double xai, yai, zai;

  Na = *na;
  Nb = *nb;

  for (i=1; i < Na; i++) 
    {
      xai = xa[i];
      yai = ya[i];
      zai = za[i];
      match[i] = 0;
      for (j=0; j < Nb; j++) 
	if(xai == xb[j] && yai == yb[j] && zai == zb[i]) {
	  match[i] = j;
	  break;
	}
    }
}

