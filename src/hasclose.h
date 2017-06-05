/*
  hasclose.h

  Function definitions to be #included in hasclose.c
  several times with different values of macros.

  Macros used:

  CLOSEFUN   name of function for pairs in a single pattern

  CROSSFUN   name of function for pairs between two patterns

  ZCOORD     if defined, coordinates are 3-dimensional

  TORUS      if defined, distances are periodic

  BUG        debugger flag

  $Revision: 1.10 $ $Date: 2017/06/05 10:53:59 $

*/

void CLOSEFUN(n,
	      x,
	      y,
#ifdef ZCOORD
	      z,
#endif
	      r,  /* distance deemed 'close' */
#ifdef TORUS
	      b,  /* box dimensions */
#endif
	      t)  /* result: true/false */
     int *n, *t;
     double *x, *y, *r;
#ifdef ZCOORD
     double *z;
#endif
#ifdef TORUS
     double *b;
#endif
{
  double xi, yi, rmax, r2max, rmaxplus, dx, dy, d2minr2;
#ifdef ZCOORD
  double zi, dz;
#endif
  int N, maxchunk, i, j;

#ifdef TORUS
  double Bx, By, Hy;
#ifdef ZCOORD
  double Bz, Hz;
#endif
#endif

  N = *n;
  rmax = *r;

  r2max = rmax * rmax;
  rmaxplus = rmax + rmax/16.0;

#ifdef TORUS
  Bx = b[0];
  By = b[1];
  Hy = By/2.0;
#ifdef ZCOORD
  Bz = b[2];
  Hz = Bz/2.0;
#endif
#endif

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
	/* scan backward from i */
	for(j = i - 1; j >= 0; j--) {
	  dx = xi - x[j];
	  if(dx > rmaxplus) 
	    break;
	  dy = y[j] - yi;
#ifdef TORUS
	  if(dy < 0.0) dy = -dy;
	  if(dy > Hy) dy = By - dy;
#endif
	  d2minr2 = dx * dx + dy * dy - r2max;
#ifdef ZCOORD
	  if(d2minr2 <= 0.0) {
	    dz = z[j] - zi;
#ifdef TORUS
	    if(dz < 0.0) dz = -dz;
	    if(dz > Hz) dz = Bz - dz;
#endif
	    d2minr2 = d2minr2 + dz * dz;
#endif
	    if(d2minr2 <= 0.0) {
	      /* pair (i, j) is close */
	      t[i] = t[j] = 1;
	    }
#ifdef ZCOORD
	  }
#endif
	}
#ifdef TORUS
	/* wrap-around */
	/* scan forward from 0 */
	for(j = 0; j < i; j++) {
	  dx = Bx + x[j] - xi;
	  if(dx > rmaxplus) 
	    break;
	  dy = y[j] - yi;
#ifdef TORUS
	  if(dy < 0.0) dy = -dy;
	  if(dy > Hy) dy = By - dy;
#endif
	  d2minr2 = dx * dx + dy * dy - r2max;
#ifdef ZCOORD
	  if(d2minr2 <= 0.0) {
	    dz = z[j] - zi;
#ifdef TORUS
	    if(dz < 0.0) dz = -dz;
	    if(dz > Hz) dz = Bz - dz;
#endif
	    d2minr2 = d2minr2 + dz * dz;
#endif
	    if(d2minr2 <= 0.0) {
	      /* pair (i, j) is close */
	      t[i] = t[j] = 1;
	    }
#ifdef ZCOORD
	  }
#endif
	}
#endif
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
#ifdef TORUS
	      b,  /* box dimensions (same for both patterns!!) */
#endif
	      t)
     int *n1, *n2, *t;
     double *x1, *y1, *x2, *y2, *r;
#ifdef ZCOORD
     double *z1, *z2;
#endif
#ifdef TORUS
     double *b;
#endif
{
  /* lengths */
  int N1, N2, maxchunk;
  /* distance parameter */
  double rmax, r2max, rmaxplus;
  /* indices */
  int i, j, jleft;
  /* temporary values */
  double x1i, y1i, xleft, dx, dy, dx2, d2minr2;
#ifdef ZCOORD
  double z1i, dz;
#endif

#ifdef TORUS
  double Bx, By, Hx, Hy;
  int jright;
#ifdef ZCOORD
  double Bz, Hz;
#endif
#endif

  N1 = *n1;
  N2 = *n2;
  rmax = *r;
  
  r2max = rmax * rmax;
  rmaxplus = rmax + rmax/16.0;

#ifdef TORUS
  Bx = b[0];
  By = b[1];
  Hx = Bx/2.0;
  Hy = By/2.0;
#ifdef BUG
  Rprintf("=> PERIODIC:  Bx = %lf, By = %lf  <= \n", Bx, By);
#endif
#ifdef ZCOORD
  Bz = b[2];
  Hz = Bz/2.0;
#endif
#endif

  if(N1 > 0 && N2 > 0) {

    i = 0; maxchunk = 0; jleft = 0;

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

#ifdef BUG
	Rprintf("------ i = %d --------\n", i);
	Rprintf(" [%d] = (%lf, %lf)\n", i, x1i, y1i);
#endif
	/* 
	   adjust starting point jleft
	*/
	xleft = x1i - rmaxplus;
	while((x2[jleft] < xleft) && (jleft+1 < N2))
	  ++jleft;

#ifdef BUG
	Rprintf("\t jleft = %d\n", jleft);
#endif

	/* 
	   process from j = jleft until dx > rmax + epsilon
	*/
	for(j=jleft; j < N2; j++) {
	  dx = x2[j] - x1i;
#ifdef BUG
	  Rprintf("\t Central loop, j = %d, dx = %lf\n", j, dx);
#endif
	  if(dx > rmaxplus)
	    break;
	  dx2 = dx * dx;
	  dy = y2[j] - y1i;
#ifdef BUG
	  Rprintf("\t\t Did not break\n\t\t dy = %lf\n", dy);
#endif
#ifdef TORUS
	  if(dy < 0.0) dy = -dy;
	  if(dy > Hy) dy = By - dy;
#ifdef BUG
	  Rprintf("\t\t periodic dy = %lf\n", dy);
#endif
#endif
	  d2minr2 = dx2 + dy * dy - r2max;
#ifdef ZCOORD
	    if(d2minr2 <= 0.0) {
	      dz = z2[j] - z1i;
#ifdef TORUS
	      if(dz < 0.0) dz = -dz;
	      if(dz > Hz) dz = Bz - dz;
#endif
	      d2minr2 = d2minr2 + dz * dz;
#endif
	      if(d2minr2 <= 0.0) {
#ifdef BUG
		Rprintf("\t\t Point %d has close neighbour\n", i);
#endif
		/* point i has a close neighbour */
		t[i] = 1;
		break;
	      }
#ifdef ZCOORD
	    }
#endif
	}

#ifdef TORUS
	jright = j;
	/* wrap-around at start */
#ifdef BUG
	Rprintf("\t Wrap around at start for j = 0 to %d\n", jleft);
#endif
	for(j=0; j < jleft; j++) {
	  dx = x1i - x2[j];
#ifdef BUG
	  Rprintf("\t\t j = %d, dx = %lf\n", j, dx);
#endif
	  if(dx < 0.0) dx = -dx;
	  if(dx > Hx) dx = Bx - dx;
#ifdef BUG
	  Rprintf("\t\t periodic dx = %lf\n", dx);
#endif
	  if(dx > rmaxplus)
	    break;
	  dx2 = dx * dx;
	  dy = y2[j] - y1i;
#ifdef BUG
	  Rprintf("\t\t Did not break\n\t\t dy = %lf\n", dy);
#endif
	  if(dy < 0.0) dy = -dy;
	  if(dy > Hy) dy = By - dy;
#ifdef BUG
	  Rprintf("\t\t periodic dy = %lf\n", dy);
#endif
	  d2minr2 = dx2 + dy * dy - r2max;
#ifdef ZCOORD
	    if(d2minr2 <= 0.0) {
	      dz = z2[j] - z1i;
	      if(dz < 0.0) dz = -dz;
	      if(dz > Hz) dz = Bz - dz;
	      d2minr2 = d2minr2 + dz * dz;
#endif
	      if(d2minr2 <= 0.0) {
		/* point i has a close neighbour */
#ifdef BUG
		Rprintf("\t\t Point %d has close neighbour\n", i);
#endif
		t[i] = 1;
		break;
	      }
#ifdef ZCOORD
	    }
#endif
	}
	/* wrap around at end */
#ifdef BUG
	Rprintf("\t Wrap around at end for j = %d to %d\n", N2-1, jright);
#endif
	for(j=N2-1; j >= jright; j--) {
	  dx = x1i - x2[j];
#ifdef BUG
	  Rprintf("\t\t j = %d, dx = %lf\n", j, dx);
#endif
	  if(dx < 0.0) dx = -dx;
	  if(dx > Hx) dx = Bx - dx;
#ifdef BUG
	  Rprintf("\t\t periodic dx = %lf\n", dx);
#endif
	  if(dx > rmaxplus)
	    break;
	  dx2 = dx * dx;
	  dy = y2[j] - y1i;
#ifdef BUG
	  Rprintf("\t\t Did not break\n\t\t dy = %lf\n", dy);
#endif
	  if(dy < 0.0) dy = -dy;
	  if(dy > Hy) dy = By - dy;
#ifdef BUG
	  Rprintf("\t\t periodic dy = %lf\n", dy);
#endif
	  d2minr2 = dx2 + dy * dy - r2max;
#ifdef ZCOORD
	    if(d2minr2 <= 0.0) {
	      dz = z2[j] - z1i;
	      if(dz < 0.0) dz = -dz;
	      if(dz > Hz) dz = Bz - dz;
	      d2minr2 = d2minr2 + dz * dz;
#endif
	      if(d2minr2 <= 0.0) {
#ifdef BUG
		Rprintf("\t\t Point %d has close neighbour\n", i);
#endif
		/* point i has a close neighbour */
		t[i] = 1;
		break;
	      }
#ifdef ZCOORD
	    }
#endif
	}
#endif
      }
    }
  }
}


