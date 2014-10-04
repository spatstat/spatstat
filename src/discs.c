#include <R.h>
#include <Rmath.h>

/*
  discs.c

  Fill binary mask with discs with given centres and radii

  $Revision: 1.3 $  $Date: 2014/10/03 10:04:36 $

*/

void discs2grid(nx, x0, xstep,  
		ny, y0, ystep,   /* pixel grid dimensions */
		nd, xd, yd, rd,  /* disc parameters */
		out)
     /* inputs */
     int *nx, *ny, *nd;
     double *x0, *xstep, *y0, *ystep;
     double *xd, *yd, *rd;
     /* output */
     int *out;
{ 
  int Nxcol, Nyrow, Ndiscs;
  int i, j, k;
  double  X0, Y0, Xstep, Ystep;
  double xk, yk, rk, rk2;
  double xj, dx, dymax; 
  int imin, imax, jmin, jmax, iminj, imaxj;

  Nxcol   = *nx;
  Nyrow   = *ny;
  Ndiscs  = *nd;
  X0      = *x0;
  Y0      = *y0;
  Xstep   = *xstep;
  Ystep   = *ystep;

  if(Ndiscs == 0)
    return;

  /* loop over discs */
  for(k = 0; k < Ndiscs; k++) {
    
    R_CheckUserInterrupt();

    xk = xd[k];
    yk = yd[k];
    rk = rd[k];

    /* find valid range of i and j */

    imax = ceil( (yk + rk - Y0)/Ystep);
    imin = floor((yk - rk - Y0)/Ystep);
    jmax = ceil( (xk + rk - X0)/Xstep);
    jmin = floor((xk - rk - X0)/Xstep);

    if(imax >= 0 && imin < Nyrow && jmax >= 0 && jmin < Nxcol) {
      
      if(imin < 0) imin = 0; 
      if(imax > Nyrow) imax = Nyrow;
      if(jmin < 0) jmin = 0; 
      if(jmax > Nxcol) jmax = Nxcol;

      rk2 = rk * rk;
      
      /* loop over relevant pixels */
      for(j = jmin; j <= jmax; j++) {

	xj = X0 + j * Xstep;
	dx = xj - xk;
	dymax = sqrt(rk2 - dx * dx);
	
	imaxj = ceil( (yk + dymax - Y0)/Ystep);
	iminj = floor((yk - dymax - Y0)/Ystep);

	if(imaxj >= 0 && iminj < Nyrow) {
	  if(iminj < 0) iminj = 0; 
	  if(imaxj > Nyrow) imaxj = Nyrow;
	  
	  for(i = iminj; i <= imaxj; i++) 
	    out[i + j * Nyrow] = 1;
	    
	}
      }
    }
  }
}





