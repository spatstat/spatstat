/*

  rasterfilter.c

  Apply linear filter to a raster image

  Copyright (C) Adrian Baddeley, Rolf Turner and Ege Rubak 2017
  Licence: GPL >= 2

  $Revision: 1.6 $  $Date: 2017/11/18 05:14:53 $


*/

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

void raster3filter(nx, ny, a, w, b)
     int *nx, *ny; /* raster dimensions */
     double *a;  /* input image */
     double *w;  /* 3x3 filter coefficients */
     double *b;  /* output image */
{ 
  int Nxcol, Nyrow, Nx1, Ny1;
  int i, j;
  double value;
  
  Nxcol   = *nx;
  Nyrow   = *ny;
  Nx1 = Nxcol - 1;
  Ny1 = Nyrow - 1;

#define A(I,J) a[(I) + (J) * Nyrow]
#define B(I,J) b[(I) + (J) * Nyrow]
#define WEIGHT(DI,DJ) w[((DI)+1) + ((DJ)+1)*3]
#define FILTER(DI,DJ) WEIGHT(DI,DJ) * A(i+(DI), j+(DJ)) 
  
  /* loop over pixels */

  for(j = 0; j < Nxcol; j++) {

    R_CheckUserInterrupt();
    
    for(i = 0; i < Nyrow; i++) {

      value = FILTER(0,0);

      if(j > 0) value += FILTER(0,-1);
      if(j < Nx1) value += FILTER(0, 1);

      if(i > 0) {
	if(j > 0) value += FILTER(-1,-1);
        value += FILTER(-1, 0);
	if(j < Nx1) value += FILTER(-1, 1);
      }
      if(i < Ny1) {
	if(j > 0) value += FILTER(1, -1);
        value += FILTER(1, 0);
	if(j < Nx1) value += FILTER(1, 1);
      }

      B(i,j) = value;
    }
  }
}



