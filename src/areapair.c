/*

  areapair.c

  $Revision: 1.6 $     $Date: 2013/09/18 04:11:42 $

  Specialised code for the second order conditional intensity
  of the area-interaction process

*/

#include <R.h>
#include <math.h>

#include "yesno.h"

/* computes area of b(A, r) \int b(B, r) \setminus \bigcup_i b(X[i], r) */

void delta2area(xa, ya, xb, yb, 
		nother, xother, yother,
		radius, epsilon, pixcount) 
     double *xa, *ya, *xb, *yb;
     int *nother;
     double *xother, *yother;
     double *radius, *epsilon;
     int *pixcount;
{ 
  int Ni, Nj, Nk, i, j, k, count, covered;
  double xA, yA, xB, yB, r, eps, r2;
  double xmin, xmax, ymin, ymax, xi, yj;
  double dxA, dyA;
  double dxB, dyB;
  double dx, dy;
  
  Nk = *nother;

  xA = *xa;
  yA = *ya;
  xB = *xb;
  yB = *yb;
  r = *radius;
  eps = *epsilon;
  r2 = r * r;

  /* find intersection of squares centred on A and B */
  if(xA < xB) {
    xmin = xB - r;
    xmax = xA + r;
  } else {
    xmin = xA - r;
    xmax = xB + r;
  }
  if(xmin > xmax) return;
  if(yA < yB) {
    ymin = yB - r;
    ymax = yA + r;
  } else {
    ymin = yA - r;
    ymax = yB + r;
  }
  if(ymin > ymax) return;
    
  /* set up grid */
  Ni = (int) ceil((xmax - xmin)/eps) + 1;
  Nj = (int) ceil((ymax - ymin)/eps) + 1;
  
  count = 0;

  for(i = 0, xi = xmin; i < Ni; i++, xi += eps) {
    dxA = xi - xA;
    for(j = 0, yj = ymin; j < Nj; j++, yj += eps) {
      dyA = yj - yA;
      if(dxA * dxA + dyA * dyA <= r2) {
	/* grid point belongs to b(A, r) */
	dxB = xi - xB;
	dyB = yj - yB;
	if(dxB * dxB + dyB * dyB <= r2) {
	  /* grid point belongs to b(A,r) \cap b(B,r) */
	  covered = NO;
	  /* test whether it is covered by another b(X[k], r) */
	  for(k = 0; k < Nk; k++) {
	    dx = xi - xother[k];
	    dy = yj - yother[k];
	    if(dx * dx + dy * dy <= r2) {
	      covered = YES;
	      break;
	    }
	  }
	  if(!covered) {
	    ++count;
	  }
	}
      }
    }
  }
  *pixcount = count;
}


