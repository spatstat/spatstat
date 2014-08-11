/*

  corrections.c

  Edge corrections

  $Revision: 1.12 $     $Date: 2013/05/27 02:09:10 $

 */

#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"
#include "yesno.h"
#include "constants.h"

#undef DEBUG


/* This constant is defined in Rmath.h */
#define TWOPI M_2PI

#define MIN(A,B) (((A) < (B)) ? (A) : (B))

#define BETWEEN(X,X0,X1) (((X) - (X0)) * ((X) - (X1)) <= 0)

#define UNDER(X,Y,X0,Y0,X1,Y1) \
  (((Y1) - (Y0)) * ((X) - (X0)) >= ((Y) - (Y0)) * ((X1)- (X0)))

#define UNDERNEATH(X,Y,X0,Y0,X1,Y1) \
    (((X0) < (X1)) ? UNDER(X,Y,X0,Y0,X1,Y1) : UNDER(X,Y,X1,Y1,X0,Y0))

#define TESTINSIDE(X,Y,X0,Y0,X1,Y1) \
  (BETWEEN(X,X0,X1) && UNDERNEATH(X, Y, X0, Y0, X1, Y1))


void ripleybox(nx, x, y, rmat, nr, xmin, ymin, xmax, ymax,  epsilon, out)
     /* inputs */
     int *nx, *nr;  /* dimensions */
     double *x, *y; /* coordinate vectors of length nx */
     double *rmat;  /* matrix nx by nr  */
     double *xmin, *ymin, *xmax, *ymax;  /* box dimensions */
     double *epsilon; /* threshold for proximity to corner */
     /* output */
     double *out;  /* output matrix nx by nr */
{
  int i, j, n, m, ijpos, ncor, maxchunk;
  double xx, yy, x0, y0, x1, y1, dL, dR, dU, dD, aL, aU, aD, aR, rij;
  double cL, cU, cD, cR, bLU, bLD, bRU, bRD, bUL, bUR, bDL, bDR;
  double corner, extang;
  double eps;

  n  = *nx;
  m  = *nr;
  x0 = *xmin;
  y0 = *ymin;
  x1 = *xmax;
  y1 = *ymax;
  eps = *epsilon;

  OUTERCHUNKLOOP(i, n, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, n, maxchunk, 16384) {
      xx = x[i];
      yy = y[i];
      /* 
	 perpendicular distance from point to each edge of rectangle
	 L = left, R = right, D = down, U = up
      */
      dL = xx - x0;
      dR = x1 - xx;
      dD = yy - y0;
      dU = y1 - yy;

      /*
	test for corner of the rectangle
      */
#define ABS(X) (((X) >= 0) ? (X) : (-X))
#define SMALL(X) ((ABS(X) < eps) ? 1 : 0)

      ncor = SMALL(dL) + SMALL(dR) + SMALL(dD) + SMALL(dU);
      corner = (ncor >= 2) ? YES : NO;
  
      /* 
	 angle between 
	 - perpendicular to edge of rectangle
	 and 
	 - line from point to corner of rectangle

      */
      bLU = atan2(dU, dL);
      bLD = atan2(dD, dL);
      bRU = atan2(dU, dR);
      bRD = atan2(dD, dR);
      bUL = atan2(dL, dU);
      bUR = atan2(dR, dU);
      bDL = atan2(dL, dD);
      bDR = atan2(dR, dD);

      for(j = 0; j < m; j++) {
	ijpos = j * n + i;
	rij = rmat[ijpos];
#ifdef DEBUG
	Rprintf("rij = %lf\n", rij);
#endif
	/*
	  half the angle subtended by the intersection between
	  the circle of radius r[i,j] centred on point i
	  and each edge of the rectangle (prolonged to an infinite line)
	*/
	aL = (dL < rij) ? acos(dL/rij) : 0.0;
	aR = (dR < rij) ? acos(dR/rij) : 0.0;
	aD = (dD < rij) ? acos(dD/rij) : 0.0;
	aU = (dU < rij) ? acos(dU/rij) : 0.0;
#ifdef DEBUG
	Rprintf("aL = %lf\n", aL);
	Rprintf("aR = %lf\n", aR);
	Rprintf("aD = %lf\n", aD);
	Rprintf("aU = %lf\n", aU);
#endif
	/* apply maxima */

	cL = MIN(aL, bLU) + MIN(aL, bLD);
	cR = MIN(aR, bRU) + MIN(aR, bRD);
	cU = MIN(aU, bUL) + MIN(aU, bUR);
	cD = MIN(aD, bDL) + MIN(aD, bDR);
#ifdef DEBUG
	Rprintf("cL = %lf\n", cL);
	Rprintf("cR = %lf\n", cR);
	Rprintf("cD = %lf\n", cD);
	Rprintf("cU = %lf\n", cU);
#endif

	/* total exterior angle over 2 pi */
	extang = (cL + cR + cU + cD)/TWOPI;

	/* add pi/2 for corners */
	if(corner) 
	  extang += 1/4;

#ifdef DEBUG
	Rprintf("extang = %lf\n", extang);
#endif
	/* OK, now compute weight */
	out[ijpos] = 1 / (1 - extang);
      }
    }
  }
}


void ripleypoly(nc, xc, yc, nr, rmat, nseg, x0, y0, x1, y1, out) 
     /* inputs */
     int *nc, *nr, *nseg;
     double *xc, *yc, *rmat;
     double *x0, *y0, *x1, *y1;
     /* output */
     double *out;
{
  int n, m, i, j, k, l, nradperpt, ncut, nchanges, maxchunk;
  double xcentre, ycentre, xx0, yy0, xx1, yy1, xx01, yy01;
  double x, y, radius, radius2, dx0, dx1, dy0;
  double a, b, c, t, det, sqrtdet, tmp;
  double theta[6], delta[7], tmid[7];
  double xtest, ytest, contrib, total;

  n = *nc;
  nradperpt = *nr;
  m = *nseg;

  OUTERCHUNKLOOP(i, n, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, n, maxchunk, 16384) {
      xcentre = xc[i];
      ycentre = yc[i];
#ifdef DEBUG
      Rprintf("centre = (%lf, %lf)\n", xcentre, ycentre);
#endif

      for(j = 0; j < nradperpt; j++) {
	radius = rmat[ j * n + i];
	radius2 = radius * radius;
#ifdef DEBUG
	Rprintf("radius = %lf\n", radius);
#endif

	total = 0.0;
	for(k=0; k < m; k++) {
#ifdef DEBUG
	  Rprintf("k = %d\n", k);
#endif
	  ncut = 0;
	  xx0 = x0[k];
	  yy0 = y0[k];
	  xx1 = x1[k];
	  yy1 = y1[k];
#ifdef DEBUG
	  Rprintf("(%lf,%lf) to (%lf,%lf)\n", xx0, yy0, xx1, yy1);
#endif
	  /* intersection with left edge */
	  dx0 = xx0 - xcentre;
	  det = radius2 - dx0 * dx0;
	  if(det > 0) {
	    sqrtdet = sqrt(det);
	    y = ycentre + sqrtdet;
	    if(y < yy0) {
	      theta[ncut] = atan2(y - ycentre, dx0);
#ifdef DEBUG
	      Rprintf("cut left at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	    y = ycentre - sqrtdet;
	    if(y < yy0) {
	      theta[ncut] = atan2(y-ycentre, dx0);
#ifdef DEBUG
	      Rprintf("cut left at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	  } else if(det == 0) {
	    if(ycentre < yy0) {
	      theta[ncut] = atan2(0.0, dx0);
#ifdef DEBUG
	      Rprintf("tangent left at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	  }
	  /* intersection with right edge */
	  dx1 = xx1 - xcentre;
	  det = radius2 - dx1 * dx1;
	  if(det > 0) {
	    sqrtdet = sqrt(det);
	    y = ycentre + sqrtdet;
	    if(y < yy1) {
	      theta[ncut] = atan2(y - ycentre, dx1);
#ifdef DEBUG
	      Rprintf("cut right at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	    y = ycentre - sqrtdet;
	    if(y < yy1) {
	      theta[ncut] = atan2(y - ycentre, dx1);
#ifdef DEBUG
	      Rprintf("cut right at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	  } else if(det == 0) {
	    if(ycentre < yy1) {
	      theta[ncut] = atan2(0.0, dx1);
#ifdef DEBUG
	      Rprintf("tangent right at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	  }
	  /* intersection with top segment */
	  xx01 = xx1 - xx0;
	  yy01 = yy1 - yy0;
	  dy0  = yy0 - ycentre;
	  a = xx01 * xx01 + yy01 * yy01;
	  b = 2 * (xx01 * dx0 + yy01 * dy0);
	  c = dx0 * dx0 + dy0 * dy0 - radius2;
	  det = b * b - 4 * a * c;
	  if(det > 0) {
	    sqrtdet = sqrt(det);
	    t = (sqrtdet - b)/(2 * a);
	    if(t >= 0 && t <= 1) {
	      x = xx0 + t * xx01;
	      y = yy0 + t * yy01;
	      theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUG
	      Rprintf("hits segment: t = %lf, theta = %lf\n", 
		      t, theta[ncut]);
#endif
	      ++ncut;
	    }
	    t = (-sqrtdet - b)/(2 * a);
	    if(t >= 0 && t <= 1) {
	      x = xx0 + t * xx01;
	      y = yy0 + t * yy01;
	      theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUG
	      Rprintf("hits segment: t = %lf, theta = %lf\n", 
		      t, theta[ncut]);
#endif
	      ++ncut;
	    }
	  } else if(det == 0) {
	    t = - b/(2 * a);
	    if(t >= 0 && t <= 1) {
	      x = xx0 + t * xx01;
	      y = yy0 + t * yy01;
	      theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUG
	      Rprintf("tangent to segment: t = %lf, theta = %lf\n", 
		      t, theta[ncut]);
#endif
	      ++ncut;
	    }
	  }
	  /* for safety, force all angles to be in range [0, 2 * pi] */
	  if(ncut > 0) 
	    for(l = 0; l < ncut; l++)
	      if(theta[l] < 0) 
		theta[l] += TWOPI;

	  /* sort angles */
	  if(ncut > 1) {
	    do {
	      nchanges = 0;
	      for(l = 0; l < ncut - 1; l++) {
		if(theta[l] > theta[l+1]) {
		  /* swap */
		  ++nchanges;
		  tmp = theta[l];
		  theta[l] = theta[l+1];
		  theta[l+1] = tmp;
		}
	      }
	    } while(nchanges > 0);
	  }
#ifdef DEBUG
	  if(ncut > 0) {
	    for(l = 0; l < ncut; l++)
	      Rprintf("theta[%d] = %lf\n", l, theta[l]);
	  }
#endif
	  /* compute length of circumference inside polygon */
	  if(ncut == 0) {
	    /* entire circle is either in or out */
	    xtest = xcentre + radius;
	    ytest = ycentre;
	    if(TESTINSIDE(xtest, ytest, xx0, yy0, xx1, yy1)) 
	      contrib = TWOPI;
	    else 
	      contrib = 0.0;
	  } else {
	    /* find midpoints and lengths of pieces (adding theta = ) */
	    delta[0] = theta[0];
	    tmid[0] = theta[0]/2;
	    if(ncut > 1) {
	      for(l = 1; l < ncut; l++) {
		delta[l] = theta[l] - theta[l-1];
		tmid[l] = (theta[l] + theta[l-1])/2;
	      }
	    }
	    delta[ncut] = TWOPI - theta[ncut - 1];
	    tmid[ncut] = (TWOPI + theta[ncut-1])/2;
	    contrib = 0.0;
	    for(l = 0; l <= ncut; l++) {
#ifdef DEBUG
	      Rprintf("delta[%d] = %lf\n", l, delta[l]);
#endif
	      xtest = xcentre + radius * cos(tmid[l]);
	      ytest = ycentre + radius * sin(tmid[l]);
	      if(TESTINSIDE(xtest, ytest, xx0, yy0, xx1, yy1)) {
		contrib += delta[l];
#ifdef DEBUG 
		Rprintf("... inside\n");
	      } else {
		Rprintf("... outside\n");
#endif
	      }

	    }
	  }
	  /* multiply by sign of trapezium */
	  if(xx0  < xx1)
	    contrib *= -1;

#ifdef DEBUG
	  Rprintf("contrib = %lf\n", contrib);
#endif
	  total += contrib;
	}
	out[ j * n + i] = total;
#ifdef DEBUG
	Rprintf("total = %lf\n", total);
#endif
      }
    }
  }
}

