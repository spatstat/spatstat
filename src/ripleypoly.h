/*

  ripleypoly.h

  Ripley's edge correction for polygonal windows

  This file is #included multiple times in corrections.c
  Macros used:

  RIPLEYFUN      Name of C function
  DEBUGPOLY      #defined if debugging information should be printed.

  TESTINSIDE     defined in corrections.c
  *CHUNKLOOP     defined in chunkloop.h
  TWOPI          defined in Rmath.h

  $Revision: 1.4 $     $Date: 2019/03/29 03:56:58 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
  Licence: GNU Public Licence >= 2

 */

void RIPLEYFUN(nc, xc, yc, nr, rmat, nseg, x0, y0, x1, y1, out) 
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
#ifdef DEBUGPOLY
      Rprintf("------- centre = (%lf, %lf) ------\n", xcentre, ycentre);
#endif

      for(j = 0; j < nradperpt; j++) {
	radius = rmat[ j * n + i];
	radius2 = radius * radius;
#ifdef DEBUGPOLY
	Rprintf("radius = %lf\n", radius);
#endif

	total = 0.0;
	for(k=0; k < m; k++) {
#ifdef DEBUGPOLY
	  Rprintf("... k = %d ... \n", k);
#endif
	  ncut = 0;
	  xx0 = x0[k];
	  yy0 = y0[k];
	  xx1 = x1[k];
	  yy1 = y1[k];
#ifdef DEBUGPOLY
	  Rprintf("(%lf,%lf) to (%lf,%lf)\n", xx0, yy0, xx1, yy1);
#endif
	  /* intersection with left edge */
	  dx0 = xx0 - xcentre;
	  det = radius2 - dx0 * dx0;
#ifdef DEBUGPOLY
	  Rprintf("Left: det = %lf\n", det);
#endif
	  if(det > 0) {
#ifdef DEBUGPOLY
	    Rprintf("\tdet > 0\n");
#endif
	    sqrtdet = sqrt(det);
	    y = ycentre + sqrtdet;
	    if(y < yy0) {
	      theta[ncut] = atan2(y - ycentre, dx0);
#ifdef DEBUGPOLY
	      Rprintf("\tcut left at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	    y = ycentre - sqrtdet;
	    if(y < yy0) {
	      theta[ncut] = atan2(y-ycentre, dx0);
#ifdef DEBUGPOLY
	      Rprintf("\tcut left at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	  } else if(det == 0) {
#ifdef DEBUGPOLY
	    Rprintf("\tdet = 0\n");
#endif
	    if(ycentre < yy0) {
	      theta[ncut] = atan2(0.0, dx0);
#ifdef DEBUGPOLY
	      Rprintf("\ttangent left at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	  }
	  /* intersection with right edge */
	  dx1 = xx1 - xcentre;
	  det = radius2 - dx1 * dx1;
#ifdef DEBUGPOLY
	  Rprintf("Right: det = %lf\n", det);
#endif
	  if(det > 0) {
#ifdef DEBUGPOLY
	    Rprintf("\tdet > 0\n");
#endif
	    sqrtdet = sqrt(det);
	    y = ycentre + sqrtdet;
	    if(y < yy1) {
	      theta[ncut] = atan2(y - ycentre, dx1);
#ifdef DEBUGPOLY
	      Rprintf("\tcut right at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	    y = ycentre - sqrtdet;
	    if(y < yy1) {
	      theta[ncut] = atan2(y - ycentre, dx1);
#ifdef DEBUGPOLY
	      Rprintf("\tcut right at theta= %lf\n", theta[ncut]);
#endif
	      ncut++;
	    }
	  } else if(det == 0) {
#ifdef DEBUGPOLY
	    Rprintf("\tdet = 0\n");
#endif
	    if(ycentre < yy1) {
	      theta[ncut] = atan2(0.0, dx1);
#ifdef DEBUGPOLY
	      Rprintf("\ttangent right at theta= %lf\n", theta[ncut]);
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
#ifdef DEBUGPOLY
	  Rprintf("Top: det = %lf\n", det);
#endif
	  if(det > 0) {
#ifdef DEBUGPOLY
	    Rprintf("\tdet > 0\n");
#endif
	    sqrtdet = sqrt(det);
	    t = (sqrtdet - b)/(2 * a);
	    if(t >= 0 && t <= 1) {
	      x = xx0 + t * xx01;
	      y = yy0 + t * yy01;
	      theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUGPOLY
	      Rprintf("\thits + segment: t = %lf, theta = %lf\n", 
		      t, theta[ncut]);
#endif
	      ++ncut;
	    }
	    t = (-sqrtdet - b)/(2 * a);
	    if(t >= 0 && t <= 1) {
	      x = xx0 + t * xx01;
	      y = yy0 + t * yy01;
	      theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUGPOLY
	      Rprintf("\thits - segment: t = %lf, theta = %lf\n", 
		      t, theta[ncut]);
#endif
	      ++ncut;
	    }
	  } else if(det == 0) {
#ifdef DEBUGPOLY
	    Rprintf("\tdet = 0\n");
#endif
	    t = - b/(2 * a);
	    if(t >= 0 && t <= 1) {
	      x = xx0 + t * xx01;
	      y = yy0 + t * yy01;
	      theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUGPOLY
	      Rprintf("\ttangent to segment: t = %lf, theta = %lf\n", 
		      t, theta[ncut]);
#endif
	      ++ncut;
	    }
	  }
#ifdef DEBUGPOLY
	  Rprintf("Finished cutting; ncut = %d\n", ncut);
#endif
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
#ifdef DEBUGPOLY
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
#ifdef DEBUGPOLY
	      Rprintf("delta[%d] = %lf\n", l, delta[l]);
#endif
	      xtest = xcentre + radius * cos(tmid[l]);
	      ytest = ycentre + radius * sin(tmid[l]);
	      if(TESTINSIDE(xtest, ytest, xx0, yy0, xx1, yy1)) {
		contrib += delta[l];
#ifdef DEBUGPOLY 
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

#ifdef DEBUGPOLY
	  Rprintf("contrib = %lf\n", contrib);
#endif
	  total += contrib;
	}
	out[ j * n + i] = total;
#ifdef DEBUGPOLY
	Rprintf("total = %lf = %lf * (2 * pi)\n\n\n", total, total/TWOPI);
#endif
      }
    }
  }
}


