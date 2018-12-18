#include <R.h>
#include <Rmath.h>

/*
  discs.c

  Fill binary mask with discs with given centres and radii

  $Revision: 1.5 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

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
  double  X0, Y0, Xstep, Ystep;

  int i, j, k;
  double xk, yk, rk, rk2, dx, dymax; 
  int imin, imax, jmin, jmax, iminj, imaxj, Nxcol1, Nyrow1;

  Nxcol   = *nx;
  Nyrow   = *ny;
  Ndiscs  = *nd;
  X0      = *x0;
  Y0      = *y0;
  Xstep   = *xstep;
  Ystep   = *ystep;

  if(Ndiscs == 0)
    return;

  Nxcol1 = Nxcol - 1;
  Nyrow1 = Nyrow - 1;

  /* loop over discs */
  for(k = 0; k < Ndiscs; k++) {
    
    R_CheckUserInterrupt();

    xk = xd[k];
    yk = yd[k];
    rk = rd[k];

    /* find valid range of i and j */

    imax = floor( (yk + rk - Y0)/Ystep);
    imin = ceil((yk - rk - Y0)/Ystep);
    jmax = floor( (xk + rk - X0)/Xstep);
    jmin = ceil((xk - rk - X0)/Xstep);

    if(imax >= 0 && imin < Nyrow && jmax >= 0 && jmin < Nxcol &&
       imax >= imin && jmax >= jmin) {
      
      if(imin < 0) imin = 0; 
      if(imax > Nyrow1) imax = Nyrow1;
      if(jmin < 0) jmin = 0; 
      if(jmax > Nxcol1) jmax = Nxcol1;

      rk2 = rk * rk;
      
      /* loop over relevant pixels */
      for(j = jmin, dx=X0 + jmin * Xstep - xk;
	  j <= jmax; 
	  j++, dx += Xstep) {

	dymax = sqrt(rk2 - dx * dx);
	
	imaxj = floor( (yk + dymax - Y0)/Ystep);
	iminj = ceil((yk - dymax - Y0)/Ystep);

	if(imaxj >= 0 && iminj < Nyrow) {
	  if(iminj < 0) iminj = 0; 
	  if(imaxj > Nyrow1) imaxj = Nyrow1;
	  
	  for(i = iminj; i <= imaxj; i++) 
	    out[i + j * Nyrow] = 1;
	    
	}
      }
    }
  }
}





