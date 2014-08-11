/*

  areadiff.c

  Area difference function

  $Revision: 1.13 $ $Date: 2012/03/27 01:38:58 $

  A(x,r) = area of disc b(0,r) not covered by discs b(x_i,r) for x_i in x
  
  Area estimated by point-counting on a fine grid

  For use in area-interaction model and related calculations

*/

#undef DEBUG

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"

/* 
   Original version areadiff()

   1 point u

   No trimming of discs

*/

void
areadiff(rad,x,y,nn,ngrid,answer) 
     /* inputs */
     double *rad;      /* radius */
     double *x, *y;    /* coordinate vectors for point pattern */
     int    *nn;       /* length of vectors x and y */
     int    *ngrid;    /* dimensions of point-counting grid */
     /* output */
     double *answer;   /* computed area */
{
  double dx, dy, xg, yg, r, r2, a2, b2, xdif, ydif;
  int i, j, k, m, n, count, covered;
  r  = *rad;
  r2 = r * r;
  n  = *nn;
  m  = *ngrid;
  dx = dy = 2 * r / (m-1);

  count = 0;

  /* run through grid points */
  for(i = 0, xg = -r; i < m; i++, xg += dx) {
    a2 = r2 - xg *xg;
    for(j = 0, yg = -r; j < m; j++, yg += dy) {
      /* test for inside disc */
      if(yg * yg < a2) {
#ifdef DEBUG
	Rprintf("\n\n (xg,yg) = (%lf, %lf)\n", xg, yg);
#endif
	/* run through data points seeking one close to (xy, yg) */
	covered = 0; 
	if(n > 0) {
	  for(k = 0; k < n; k++) {
#ifdef DEBUG
	    Rprintf("(x[%d],y[%d]) = (%lf,%lf)\n", k, k, x[k], y[k]);
#endif
	    xdif = x[k] - xg;
	    b2 = r2 - xdif * xdif;
	    if(b2 > 0) {
	      ydif = y[k] - yg;
	      if(b2 - ydif * ydif > 0) {
#ifdef DEBUG
		Rprintf("(x[%d], y[%d]) = (%lf, %lf) covers!\n", 
			k, k, x[k], y[k]);
#endif
		covered = 1;
		break;
	      }
	    }
	  }
	}
	if(covered == 0) {
	  ++count;
#ifdef DEBUG
	  Rprintf("Not covered; incrementing count\n");
#endif
	}
      }
    }
  }

#ifdef DEBUG
  Rprintf("Count = %d\n", count);
#endif
  
  /* calculate area */
  *answer = ((double) count) * dx * dy;
}

/* similar function, handles multiple values of 'r' */

void
areadifs(rad,nrads,x,y,nxy,ngrid,answer) 
     /* inputs */
     double *rad;      /* vector of radii */
     int    *nrads;     /* length of 'rads' */
     double *x, *y;    /* coordinate vectors for point pattern */
     int    *nxy;       /* length of vectors x and y */
     int    *ngrid;    /* dimensions of point-counting grid */
     /* output */
     double *answer;   /* computed areas (vector of length 'nrads') */
{
  double dx, dy, xg, yg, xg2, r, r2, a2, b2, xdif, ydif;
  int i, j, k, l, m, n, nr, m0, count, covered, maxchunk;

  n  = *nxy;
  nr = *nrads;
  m  = *ngrid;

  /* run through radii in chunks of 2^14 */
  OUTERCHUNKLOOP(l, nr, maxchunk, 16384) {

    R_CheckUserInterrupt();

    INNERCHUNKLOOP(l, nr, maxchunk, 16384) {
      r  = rad[l];
      if(r == 0.0) {
	answer[l] = 0.0;
      } else if(n == 0) {
	answer[l] = M_PI * r * r;
      } else {
	r2 = r * r;
	dx = dy = 2 * r / (m-1);
	count = 0;

	/* run through grid points in disc of radius r */
	for(i = 0, xg = -r; i < m; i++, xg += dx) {
	  a2 = r2 - xg * xg;
	  m0 = (a2 > 0.0) ? floor(sqrt(a2)/dy) : 0;
	  for(j = -m0, yg = -m0 * dy; j <= m0; j++, yg += dy) {
#ifdef DEBUG
	    Rprintf("\n\n (xg,yg) = (%lf, %lf)\n", xg, yg);
#endif
	    /* run through data points seeking one close to (xy, yg) */
	    covered = 0;
	    for(k = 0; k < n; k++) {
#ifdef DEBUG
	      Rprintf("(x[%d],y[%d]) = (%lf,%lf)\n", k, k, x[k], y[k]);
#endif
	      xdif = x[k] - xg;
	      b2 = r2 - xdif * xdif;
	      if(b2 > 0) {
		ydif = y[k] - yg;
		if(b2 - ydif * ydif > 0) {
#ifdef DEBUG
		  Rprintf("(x[%d], y[%d]) = (%lf, %lf) covers!\n", 
			  k, k, x[k], y[k]);
#endif
		  covered = 1;
		  break;
		}
	      }
	    } /* end of loop through data points */
	    if(covered == 0) {
	      ++count;
#ifdef DEBUG
	      Rprintf("Not covered; incrementing count\n");
#endif
	    }
	  }
	} /* end of loop over grid points */

#ifdef DEBUG
	Rprintf("Count = %d\n", count);
#endif
  
	/* calculate area for this value of r*/
	answer[l] = ((double) count) * dx * dy;
      }
      /* end of if(r==0).. else {...} */
    }
  }
}

/*
    Modified version

    multiple test points u
    
    discs constrained inside a rectangle

*/

void
areaBdif(rad,nrads,x,y,nxy,ngrid,x0,y0,x1,y1,answer) 
     /* inputs */
     double *rad;      /* vector of radii */
     int    *nrads;     /* length of 'rads' */
     double *x, *y;    /* coordinate vectors for point pattern */
     int    *nxy;       /* length of vectors x and y */
     int    *ngrid;    /* dimensions of point-counting grid */
     double *x0,*y0,*x1,*y1;  /* constraint rectangle */
     /* output */
     double *answer;   /* computed areas (vector of length 'nrads') */
{
  double dx, dy, xg, yg, r, r2, a, a2, b2, xdif, ydif;
  double xleft, xright, ylow, yhigh;
  double xmin, ymin, xmax, ymax;
  int i, j, k, l, m, n, nr, ileft, iright, mlow, mhigh, count, covered;

  n  = *nxy;
  nr = *nrads;
  m  = *ngrid;

  xmin = *x0;
  ymin = *y0;
  xmax = *x1;
  ymax = *y1;

  /* run through radii */
  for(l = 0; l < nr; l++) {
    r  = rad[l];
    if(r == 0.0) {
      answer[l] = 0.0;
    } else if (n == 0) {
      answer[l]= M_PI * r * r;
    } else {
      r2 = r * r;
      dx = dy = 2 * r / (m-1);
      count = 0;

      /* run through grid points in disc intersected with box */
      xleft = (xmin > -r) ? xmin : -r;
      xright = (xmax < r) ? xmax : r;
      ileft = ceil(xleft/dx);
      iright = floor(xright/dx);

      if(ileft <= iright) {
	for(i = ileft, xg = ileft * dx; i <= iright; i++, xg += dx) {
	  a2 = r2 - xg * xg;
	  a = (a2 > 0) ? sqrt(a2): 0.0;
	  yhigh = (ymax < a) ? ymax: a;
	  ylow  = (ymin > -a) ? ymin: -a;
	  mhigh = floor(yhigh/dy);
	  mlow  = ceil(ylow/dy);
	  if(mlow <= mhigh) {
	    for(j = mlow, yg = mlow * dy; j <= mhigh; j++, yg += dy) {
#ifdef DEBUG
	      Rprintf("\n\n (xg,yg) = (%lf, %lf)\n", xg, yg);
#endif
	      /* run through data points seeking one close to (xy, yg) */
	      covered = 0;
	      for(k = 0; k < n; k++) {
#ifdef DEBUG
		Rprintf("(x[%d],y[%d]) = (%lf,%lf)\n", 
			k, k, x[k], y[k]);
#endif
		xdif = x[k] - xg;
		b2 = r2 - xdif * xdif;
		if(b2 > 0) {
		  ydif = y[k] - yg;
		  if(b2 - ydif * ydif > 0) {
#ifdef DEBUG
		    Rprintf("(x[%d], y[%d]) = (%lf, %lf) covers!\n", 
			    k, k, x[k], y[k]);
#endif
		    covered = 1;
		    break;
		  }
		}
	      }
	      /* end of loop over data points */
	      if(covered == 0) {
		++count;
#ifdef DEBUG
		Rprintf("Not covered; incrementing count\n");
#endif
	      }
	    }
	  }
	}
      }
      /* end of loop over grid points */

#ifdef DEBUG
      Rprintf("Count = %d\n", count);
#endif
  
      /* calculate area for this value of r*/
      answer[l] = ((double) count) * dx * dy;
    }
    /* end of if(r==0).. else {...} */
  }
  /* end of loop over r values */
}




