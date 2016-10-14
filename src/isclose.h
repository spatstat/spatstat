/*
  isclose.h

  Function definitions to be #included in isclose.c
  several times with different values of macros.

  Macros used:

  CLOSEFUN   name of function for pairs in a single pattern

  CROSSFUN   name of function for pairs between two patterns

  ZCOORD     if defined, coordinates are 3-dimensional

  $Revision: 1.3 $ $Date: 2016/10/14 07:41:03 $

*/

void CLOSEFUN(n,
	      x,
	      y,
#ifdef ZCOORD
	      z,
#endif
	      r,  /* distance deemed 'close' */
	      t)  /* result: true/false */
     int *n, *t;
     double *x, *y, *r;
#ifdef ZCOORD
     double *z;
#endif
{
  double xi, yi, rmax, r2max, rmaxplus, dx, dy, d2;
#ifdef ZCOORD
  double zi, dz;
#endif
  int N, maxchunk, i, j;

  N = *n;
  rmax = *r;
  
  r2max = rmax * rmax;
  rmaxplus = rmax + rmax/16.0;

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < N) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > N) maxchunk = N;

    for(; i < maxchunk; i++) {

      xi = x[i];
      yi = y[i];
#ifdef ZCOORD
      zi = z[i];
#endif

      if(i > 0) {
	/* scan backward */
	for(j = i - 1; j >= 0; j--) {
	  dx = xi - x[j];
	  if(dx > rmaxplus) 
	    break;
	  dy = y[j] - yi;
	  d2 = dx * dx + dy * dy;
#ifdef ZCOORD
	  if(d2 <= r2max) {
	    dz = z[j] - zi;
	    d2 = d2 + dz * dz;
#endif
	    if(d2 <= r2max) {
	      /* pair (i, j) is close */
	      t[i] = t[j] = 1;
	    }
#ifdef ZCOORD
	  }
#endif
	}
      }
    }
  }
}

/* ........................................................ */

void CROSSFUN(n1,
	      x1,
	      y1,
#ifdef ZCOORD
	      z1,
#endif
	      n2,
	      x2,
	      y2,
#ifdef ZCOORD
	      z2,
#endif
	      r,
	      t)
     int *n1, *n2, *t;
     double *x1, *y1, *x2, *y2, *r;
#ifdef ZCOORD
     double *z1, *z2;
#endif
{
  /* lengths */
  int N1, N2, maxchunk;
  /* distance parameter */
  double rmax, r2max, rmaxplus;
  /* indices */
  int i, j, jleft;
  /* temporary values */
  double x1i, y1i, xleft, dx, dy, dx2, d2;
#ifdef ZCOORD
  double z1i, dz;
#endif

  N1 = *n1;
  N2 = *n2;
  rmax = *r;
  
  r2max = rmax * rmax;
  rmaxplus = rmax + rmax/16.0;

  if(N1 > 0 && N2 > 0) {

    i = 0; maxchunk = 0;

    while(i < N1) {

      R_CheckUserInterrupt();

      maxchunk += 65536;
      if(maxchunk > N1) maxchunk = N1;

      for( ; i < maxchunk; i++) {

	x1i = x1[i];
	y1i = y1[i];
#ifdef ZCOORD
	z1i = z1[i];
#endif

	/* 
	   adjust starting point jleft
	*/
	xleft = x1i - rmaxplus;
	while((x2[jleft] < xleft) && (jleft+1 < N2))
	  ++jleft;

	/* 
	   process from j = jleft until dx > rmax + epsilon
	*/
	for(j=jleft; j < N2; j++) {

	  /* squared interpoint distance */
	  dx = x2[j] - x1i;
	  if(dx > rmaxplus)
	    break;
	  dx2 = dx * dx;
	  dy = y2[j] - y1i;
	  d2 = dx2 + dy * dy;
#ifdef ZCOORD
	    if(d2 <= r2max) {
	      dz = z2[j] - z1i;
	      d2 = d2 + dz * dz;
#endif
	      if(d2 <= r2max) {
		/* point i has a close neighbour */
		t[i] = 1;
		break;
	      }
#ifdef ZCOORD
	    }
#endif
	}
      }
    }
  }
}

