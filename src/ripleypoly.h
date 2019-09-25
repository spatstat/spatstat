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

  $Revision: 1.20 $     $Date: 2019/09/25 05:58:56 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
  Licence: GNU Public Licence >= 2

 */

#undef DEBUGLEVEL

#ifndef DEBUGPOLY
#define DEBUGLEVEL 0
#else
#define DEBUGLEVEL 3
#endif

/* SPLITPOINT is used only when DEBUGLEVEL = 2 */
#undef SPLITPOINT
#define SPLITPOINT 0

#undef ROUNDED
#ifdef _WIN32
/* Avoid quirks of Windows i386 */
#define ROUNDED(X) ((float)(X))
#else
#define ROUNDED(X) ((float)(X)) 
/* WAS: define ROUNDED(X) ((double)(X)) */
#endif

void RIPLEYFUN(nc, xc, yc, bd, nr, rmat, nseg, x0, y0, x1, y1, out) 
     /* inputs */
     int *nc, *nr, *nseg;
     double *xc, *yc, *bd, *rmat;
     double *x0, *y0, *x1, *y1;
     /* output */
     double *out;
{
  int n, m, i, j, k, l, nradperpt, ncut, nchanges, maxchunk;
  double xcentre, ycentre, xx0, yy0, xx1, yy1, xx01, yy01;
  double bdisti;
  double x, y, radius, radius2, dx0, dx1, dy0;
  double a, b, c, t, det, sqrtdet, tmp;
  double theta[6], delta[7], tmid[7];
  double xtest, ytest, contrib, total;

  n = *nc;
  nradperpt = *nr;
  m = *nseg;

#if (DEBUGLEVEL == 2)  
  Rprintf("/// Debug level 2, split point %d ///\n", (int) SPLITPOINT);
#elif (DEBUGLEVEL > 0)
  Rprintf("/// Debug level %d ///\n", (int) DEBUGLEVEL);
#endif
  
  OUTERCHUNKLOOP(i, n, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, n, maxchunk, 16384) {
      xcentre = xc[i];
      ycentre = yc[i];
      bdisti  = bd[i];
#if (DEBUGLEVEL >= 3)
      Rprintf("------- centre[%d] = (%lf, %lf) ------\n", i, xcentre, ycentre);
      Rprintf("        boundary distance %lf \n", bdisti);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 1))
      Rprintf("------- centre[%d] ------\n", i);
#endif      
      for(j = 0; j < nradperpt; j++) {
	radius = rmat[ j * n + i];
	radius2 = (double) (radius * radius);
#if (DEBUGLEVEL >= 3)	
	Rprintf("radius[%d] = %lf\n", j, radius);
#elif (DEBUGLEVEL >= 2)	
	Rprintf("radius[%d]\n", j);
#endif
	if(bdisti > radius) {
	  /* no crossings */
	  total = TWOPI;
#if (DEBUGLEVEL >= 2)	
	  Rprintf("no crossings; total = 2*pi\n");
#endif
	} else {
	  /* run through all boundary segments */
	  total = 0.0;
	  for(k=0; k < m; k++) {
	    ncut = 0;
	    xx0 = x0[k];
	    yy0 = y0[k];
	    xx1 = x1[k];
	    yy1 = y1[k];
#if (DEBUGLEVEL >= 3)	  
	    Rprintf("... Edge[%d] = (%lf,%lf) to (%lf,%lf)\n",
		    k, xx0, yy0, xx1, yy1);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 2))
	    Rprintf("... Edge[%d]\n", k);
#endif
	    /* intersection with left edge */
	    dx0 = xx0 - xcentre;
	    det = (double) (radius2 - dx0 * dx0);
#if (DEBUGLEVEL >= 3)	  
	    Rprintf("Left: det = %lf\n", det);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 3))
	    Rprintf("Left:\n");
#endif
	    if(ROUNDED(det) > ROUNDED(0.0)) {
#if ((DEBUGLEVEL >= 3) || ((DEBUGLEVEL == 2) && (SPLITPOINT >= 4)))
	      Rprintf("\tdet > 0\n");
#endif	    
	      sqrtdet = (double) sqrt(det);
	      y = (double) (ycentre + sqrtdet);
	      if(ROUNDED(y) < ROUNDED(yy0)) {
		theta[ncut] = (double) atan2(y - ycentre, dx0);
#if (DEBUGLEVEL >= 3)	      
		Rprintf("\tcut left at theta= %lf\n", theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 5))
		Rprintf("\tcut left (+)\n");
#endif
		ncut++;
	      }
	      y = (double) (ycentre - sqrtdet);
	      if(ROUNDED(y) < ROUNDED(yy0)) {
		theta[ncut] = (double) atan2(y-ycentre, dx0);
#if (DEBUGLEVEL >= 3)	      
		Rprintf("\tcut left at theta= %lf\n", theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 6))
		Rprintf("\tcut left (-)\n");
#endif
		ncut++;
	      }
	    } else if(ROUNDED(det) == ROUNDED(0.0)) {
#if ((DEBUGLEVEL >= 3) || ((DEBUGLEVEL == 2) && (SPLITPOINT >= 7)))
	      Rprintf("\tdet = 0\n");
#endif
	      if(ROUNDED(ycentre) < ROUNDED(yy0)) {
		theta[ncut] = (double) atan2(0.0, dx0);
#if (DEBUGLEVEL >= 3)	      
		Rprintf("\ttangent left at theta= %lf\n", theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 8))
		Rprintf("\ttangent left\n");
#endif
		ncut++;
	      }
	    }
	    /* intersection with right edge */
	    dx1 = xx1 - xcentre;
	    det = (double) (radius2 - dx1 * dx1);
#if (DEBUGLEVEL >= 3)	  
	    Rprintf("Right: det = %lf\n", det);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 9))
	    Rprintf("Right:\n");
#endif
	    if(ROUNDED(det) > ROUNDED(0.0)) {
#if ((DEBUGLEVEL >= 3) || ((DEBUGLEVEL == 2) && (SPLITPOINT >= 10)))
	      Rprintf("\tdet > 0\n");
#endif
	      sqrtdet = (double) sqrt(det);
	      y = (double) (ycentre + sqrtdet);
	      if(ROUNDED(y) < ROUNDED(yy1)) {
		theta[ncut] = (double) atan2(y - ycentre, dx1);
#if (DEBUGLEVEL >= 3)	      
		Rprintf("\tcut right at theta= %lf\n", theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 11))
		Rprintf("\tcut right (+)\n");
#endif
		ncut++;
	      }
	      y = (double) (ycentre - sqrtdet);
	      if(ROUNDED(y) < ROUNDED(yy1)) {
		theta[ncut] = (double) atan2(y - ycentre, dx1);
#if (DEBUGLEVEL >= 3)	      
		Rprintf("\tcut right at theta= %lf\n", theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 12))
		Rprintf("\tcut right (-)\n");
#endif
		ncut++;
	      }
	    } else if(ROUNDED(det) == ROUNDED(0.0)) {
#if ((DEBUGLEVEL >= 3) || ((DEBUGLEVEL == 2) && (SPLITPOINT >= 13)))
	      Rprintf("\tdet = 0\n");
#endif
	      if(ycentre < yy1) {
		theta[ncut] = (double) atan2(0.0, dx1);
#if (DEBUGLEVEL >= 3)	      
		Rprintf("\ttangent right at theta= %lf\n", theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 14))
		Rprintf("\ttangent right\n");
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
	    det = (double) (b * b - 4 * a * c);
#if (DEBUGLEVEL >= 3)	  
	    Rprintf("Top: det = %lf\n", det);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 15))
	    Rprintf("Top:\n");
#endif	  
	    if(ROUNDED(det) > ROUNDED(0.0)) {
#if ((DEBUGLEVEL >= 3) || ((DEBUGLEVEL == 2) && (SPLITPOINT >= 16)))
	      Rprintf("\tdet > 0\n");
#endif
	      sqrtdet = (double) sqrt(det);
	      t = (double) ((sqrtdet - b)/(2 * a));
	      if(ROUNDED(0.0) <= ROUNDED(t) && ROUNDED(t) <= ROUNDED(1.0)) {
		x = xx0 + t * xx01;
		y = yy0 + t * yy01;
		theta[ncut] = (double) atan2(y - ycentre, x - xcentre);
#if (DEBUGLEVEL >= 3)	      
		Rprintf("\thits + segment: t = %lf, theta = %lf\n", 
			t, theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 17))
		Rprintf("\thits + segment\n");
#endif
		++ncut;
	      }
	      t = (double) ((-sqrtdet - b)/(2 * a));
	      if(ROUNDED(0.0) <= ROUNDED(t) && ROUNDED(t) <= ROUNDED(1.0)) {
		x = xx0 + t * xx01;
		y = yy0 + t * yy01;
		theta[ncut] = (double) atan2(y - ycentre, x - xcentre);
#if (DEBUGLEVEL >= 3)	      	      
		Rprintf("\thits - segment: t = %lf, theta = %lf\n", 
			t, theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 18))
		Rprintf("\thits - segment\n");
#endif
		++ncut;
	      }
	    } else if(ROUNDED(det) == ROUNDED(0.0)) {
#if ((DEBUGLEVEL >= 3) || ((DEBUGLEVEL == 2) && (SPLITPOINT >= 19)))
	      Rprintf("\tdet = 0\n");
#endif
	      t = (double) (- b/(2 * a));
	      if(ROUNDED(0.0) <= ROUNDED(t) && ROUNDED(t) <= ROUNDED(1.0)) {
		x = xx0 + t * xx01;
		y = yy0 + t * yy01;
		theta[ncut] = (double) atan2(y - ycentre, x - xcentre);
#if (DEBUGLEVEL >= 3)	      	      	      
		Rprintf("\ttangent to segment: t = %lf, theta = %lf\n", 
			t, theta[ncut]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 20))
		Rprintf("\ttangent to segment\n");
#endif
		++ncut;
	      }
	    }
#if ((DEBUGLEVEL >= 3) || ((DEBUGLEVEL == 2) && (SPLITPOINT >= 21)))
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
#if (DEBUGLEVEL >= 3)	  
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
#if (DEBUGLEVEL >= 3)	      
		Rprintf("Interval %d, width %lf:", l, delta[l]);
#elif ((DEBUGLEVEL == 2) && (SPLITPOINT >= 22))
		Rprintf("Interval %d:", l);
#endif
		xtest = (double) (xcentre + radius * cos(tmid[l]));
		ytest = (double) (ycentre + radius * sin(tmid[l]));
		if(TESTINSIDE(xtest, ytest, xx0, yy0, xx1, yy1)) {
		  contrib += delta[l];
#if ((DEBUGLEVEL >= 3) || ((DEBUGLEVEL == 2) && (SPLITPOINT >= 23)))
		  Rprintf("inside\n");
		} else {
		  Rprintf("outside\n");
#endif
		}
	      }
	    }
	    /* multiply by sign of trapezium */
	    if(xx0  < xx1)
	      contrib = -contrib;

#if (DEBUGLEVEL >= 3)
	    Rprintf("contrib = %lf\n", contrib);
#endif
	    total += contrib;
	  }
	}
	out[ j * n + i] = total;
#if (DEBUGLEVEL >= 1)
	Rprintf("\nTotal = %lf = %lf * (2 * pi)\n", total, total/TWOPI);
#endif
      }
    }
  }
}
