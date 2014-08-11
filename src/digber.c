/*

  digber.c

  Diggle-Berman function J used in bandwidth selection

  J(r) = \int_0^(2r) phi(t, r) dK(t)

  where K is the K-function and
       phi(t, r) = 2 r^2 * (acos(y) - y sqrt(1 - y^2))
                        where y = t/(2r).

  $Revision: 1.7 $     $Date: 2013/08/24 11:13:43 $

 */

#include <math.h>

double sqrt(), acos();

/* 
   r is the vector of distance values, starting from 0, with length nr,
   equally spaced.

   dK = diff(K) is the vector of increments of the K-function,
   with length ndK = nr-1.

   values of J are computed only up to max(r)/2

   nrmax = floor(nr/2).

*/

void digberJ(r, dK, nr, nrmax, ndK, J) 
     /* inputs */
     int *nr, *nrmax, *ndK;
     double *r, *dK;
     /* output */
     double *J;  
{ 
  int i, j, Ni, NdK;
  double ri, twori, tj, y, phiy, integral;

  Ni = *nrmax;
  NdK = *ndK;

  J[0] = 0.0;

  for(i = 1; i < Ni; i++) {
    ri = r[i];
    twori = 2 * ri;
    integral = 0.0;
    for(j = 0; j < NdK; j++) {
      tj = r[j];
      y = tj/twori;
      if(y >= 1.0) break;
      phiy = acos(y) - y * sqrt(1 - y * y);
      integral += phiy * dK[j];
    }
    J[i] = 2 * ri * ri * integral;
  }
}


  

  
