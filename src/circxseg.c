/*

  circxseg.c

  Intersections between circles and line segments

  circXseg: centres * radii * segments 

  circMseg: matrix of radii with rows corresponding to centres

  circPseg: parallel vectors of centres and radii

  $Revision: 1.8 $     $Date: 2017/01/21 10:50:32 $

 */

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

#undef BUGGER

/* circXseg: consider every combination of centre, radius, segment */

SEXP circXseg(SEXP XC,        /* circle centres */
              SEXP YC,        
	      SEXP R,         /* radii */
	      SEXP X0,        /* segments */
              SEXP Y0, 
              SEXP X1, 
	      SEXP Y1         
	      )
{
  double *xc, *yc, *r, *x0, *y0, *x1, *y1;
  int Nc, Ns, Nr, Ne, NeMax, newmax;

  /* outputs */
  int *ie, *je, *ke;  /* provenance of each intersection */
  double *xe, *ye, *sinalpha; /* cut points and angles */
  SEXP out, iout, jout, kout, xout, yout, sinout;
  int *ip, *jp, *kp;
  double *xp, *yp, *sp;

  /* internal */
  int i, j, k, m;
  double xci, yci, rk, x0c, y0c, dx, dy, A, B, C, Det, sqrtDet, sina;
  double u, u1, u2, slope, intcept, xcut, ycut, xnorm, ynorm, hx, hy;
  double twoA, fourA, Bsquared, Cdist2;

  PROTECT(XC = AS_NUMERIC(XC));
  PROTECT(YC = AS_NUMERIC(YC));
  PROTECT(R  = AS_NUMERIC(R));
  PROTECT(X0 = AS_NUMERIC(X0));
  PROTECT(Y0 = AS_NUMERIC(Y0));
  PROTECT(X1 = AS_NUMERIC(X1));
  PROTECT(Y1 = AS_NUMERIC(Y1));
  /* That's 7 protected */

  /* get pointers */
  xc = NUMERIC_POINTER(XC);
  yc = NUMERIC_POINTER(YC);
  r  = NUMERIC_POINTER(R);
  x0 = NUMERIC_POINTER(X0);
  y0 = NUMERIC_POINTER(Y0);
  x1 = NUMERIC_POINTER(X1);
  y1 = NUMERIC_POINTER(Y1);

  /* determine lengths of vectors */
  Nc = LENGTH(XC);
  Nr = LENGTH(R);
  Ns = LENGTH(X0);

  /* Guess amount of storage required */
  NeMax = 4 * (Ns + Nr + Nc);
  if(NeMax < 2048) NeMax = 2048;
  ie = (int *) R_alloc(NeMax, sizeof(int));
  je = (int *) R_alloc(NeMax, sizeof(int));
  ke = (int *) R_alloc(NeMax, sizeof(int));
  xe = (double *) R_alloc(NeMax, sizeof(double));
  ye = (double *) R_alloc(NeMax, sizeof(double));
  sinalpha = (double *) R_alloc(NeMax, sizeof(double));

  /* initialise output */
  Ne = 0;

  if(Nc > 0 && Ns > 0 && Nr > 0) {
    /* loop over circle centres */
    for(i = 0; i < Nc; i++) {
#ifdef BUGGER
      Rprintf("Circle %d\n", i);
#endif
      R_CheckUserInterrupt();
      xci = xc[i];
      yci = yc[i];
      /* loop over segments */
      for(j = 0; j < Ns; j++) {
#ifdef BUGGER
	Rprintf("\tSegment %d\n", j);
#endif
	dx = x1[j] - x0[j];
	dy = y1[j] - y0[j];
	x0c = x0[j] - xci;
	y0c = y0[j] - yci;
	/* find intersections between circle and infinite line */
	A = dx * dx + dy * dy;
	B = 2 * (dx * x0c + dy * y0c);
	twoA = 2.0 * A;
	fourA = 4.0 * A;
	Bsquared = B * B;
	Cdist2 = x0c * x0c + y0c * y0c;
	/* loop over radii */
	for(k = 0; k < Nr; k++) {
#ifdef BUGGER
	  Rprintf("\t\tRadius %d\n", k);
#endif
	  rk = r[k];
	  C = Cdist2 - rk * rk;
	  Det = Bsquared - fourA * C;
	  if(Det > 0.0) {
	    /* two intersection points */
	    sqrtDet = sqrt(Det);
	    u1 = (-B - sqrtDet)/twoA;
	    u2 = (-B + sqrtDet)/twoA;
	    if(u1 > 0.0 && u1 < 1.0) {
	      /* first intersection point is inside segment */
	      if(dx != 0) {
		/* sloping line */
		slope = dy/dx;
		intcept = y0c - slope * x0c;
		xcut = x0c + u1 * dx;
		ycut = y0c + u1 * dy;
		ynorm = intcept/(slope * slope + 1.0);
		xnorm = - ynorm * slope;
	      } else {
		/* vertical line */
		xcut = x0c;
		ycut = y0c + u1 * dy;
		xnorm = x0c;
		ynorm = 0.0;
	      }
	      hx = xcut - xnorm;
	      hy = ycut - ynorm;

	      sina = sqrt(hx * hx + hy * hy)/rk;
	      if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	      /* add to output */
#ifdef BUGGER
	      Rprintf("\t\t\tAdding..\n");
#endif
	      sinalpha[Ne] = sina;
	      xe[Ne] = xcut + xci;
	      ye[Ne] = ycut + yci;
	      ie[Ne] = i + 1;
	      je[Ne] = j + 1;
	      ke[Ne] = k + 1;
	      ++Ne;
	      if(Ne >= NeMax) {
		/* storage overflow; reallocate */
#ifdef BUGGER
		Rprintf("\t\t\tOVERFLOW..\n");
#endif
		newmax = 2 * NeMax;
		xe = (double *) S_realloc((char *) xe,
					  newmax, NeMax, sizeof(double));
		ye = (double *) S_realloc((char *) ye,
					  newmax, NeMax, sizeof(double));
		ie = (int *) S_realloc((char *) ie,
				       newmax, NeMax, sizeof(int));
		je = (int *) S_realloc((char *) je,
				       newmax, NeMax, sizeof(int));
		ke = (int *) S_realloc((char *) ke,
				       newmax, NeMax, sizeof(int));
		sinalpha = (double *) S_realloc((char *) sinalpha,
						newmax, NeMax, sizeof(double));
		NeMax = newmax;
	      }
	    }
	    if(u2 > 0.0 && u2 < 1.0) {
	      /* second intersection point is inside segment */
	      if(dx != 0) {
		/* sloping line */
		slope = dy/dx;
		intcept = y0c - slope * x0c;
		xcut = x0c + u2 * dx;
		ycut = y0c + u2 * dy;
		ynorm = intcept/(slope * slope + 1.0);
		xnorm = - ynorm * slope;
	      } else {
		/* vertical line */
		xcut = x0c;
		ycut = y0c + u2 * dy;
		xnorm = x0c;
		ynorm = 0.0;
	      }
	      hx = xcut - xnorm;
	      hy = ycut - ynorm;

	      sina = sqrt(hx * hx + hy * hy)/rk;
	      if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	      /* add to output */
#ifdef BUGGER
	      Rprintf("\t\t\tAdding..\n");
#endif
	      sinalpha[Ne] = sina;
	      xe[Ne] = xcut + xci;
	      ye[Ne] = ycut + yci;
	      ie[Ne] = i + 1;
	      je[Ne] = j + 1;
	      ke[Ne] = k + 1;
	      ++Ne;
	      if(Ne >= NeMax) {
		/* storage overflow; reallocate */
#ifdef BUGGER
		Rprintf("\t\t\tOVERFLOW..\n");
#endif
		newmax = 2 * NeMax;
		xe = (double *) S_realloc((char *) xe,
					  newmax, NeMax, sizeof(double));
		ye = (double *) S_realloc((char *) ye,
					  newmax, NeMax, sizeof(double));
		ie = (int *) S_realloc((char *) ie,
				       newmax, NeMax, sizeof(int));
		je = (int *) S_realloc((char *) je,
				       newmax, NeMax, sizeof(int));
		ke = (int *) S_realloc((char *) ke,
				       newmax, NeMax, sizeof(int));
		sinalpha = (double *) S_realloc((char *) sinalpha,
						newmax, NeMax, sizeof(double));
		NeMax = newmax;
	      }
	    }
	  } else if(Det == 0.0) {
	    /* tangent point */
	    u = -B/twoA;
	    if(u > 0.0 && u < 1.0) {
	      /* tangent point is inside segment */
	      if(dx != 0) {
		/* sloping line */
		slope = dy/dx;
		intcept = y0c - slope * x0c;
		xcut = x0c + u * dx;
		ycut = y0c + u * dy;
		ynorm = intcept/(slope * slope + 1.0);
		xnorm = - ynorm * slope;
	      } else {
		/* vertical line */
		xcut = x0c;
		ycut = y0c + u * dy;
		xnorm = x0c;
		ynorm = 0.0;
	      }
	      hx = xcut - xnorm;
	      hy = ycut - ynorm;

	      sina = sqrt(hx * hx + hy * hy)/rk;
	      if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	      /* add to output */
#ifdef BUGGER
	      Rprintf("\t\t\tAdding..\n");
#endif
	      sinalpha[Ne] = sina;
	      xe[Ne] = xcut + xci;
	      ye[Ne] = ycut + yci;
	      ie[Ne] = i + 1;
	      je[Ne] = j + 1;
	      ke[Ne] = k + 1;
	      ++Ne;
	      if(Ne >= NeMax) {
		/* storage overflow; reallocate */
#ifdef BUGGER
		Rprintf("\t\t\tOVERFLOW..\n");
#endif
		newmax = 2 * NeMax;
		xe = (double *) S_realloc((char *) xe,
					  newmax, NeMax, sizeof(double));
		ye = (double *) S_realloc((char *) ye,
					  newmax, NeMax, sizeof(double));
		ie = (int *) S_realloc((char *) ie,
				       newmax, NeMax, sizeof(int));
		je = (int *) S_realloc((char *) je,
				       newmax, NeMax, sizeof(int));
		ke = (int *) S_realloc((char *) ke,
				       newmax, NeMax, sizeof(int));
		sinalpha = (double *) S_realloc((char *) sinalpha,
						newmax, NeMax, sizeof(double));
		NeMax = newmax;
	      }
	    }
	  }
	}
      }
    }
  }

  /* pack up */
  PROTECT(out = NEW_LIST(6));
  PROTECT(iout = NEW_INTEGER(Ne));
  PROTECT(jout = NEW_INTEGER(Ne));
  PROTECT(kout = NEW_INTEGER(Ne));
  PROTECT(xout = NEW_NUMERIC(Ne));  
  PROTECT(yout = NEW_NUMERIC(Ne));  
  PROTECT(sinout = NEW_NUMERIC(Ne));  
  /* 7 + 1 + 6 = 14 protected */
  ip = INTEGER_POINTER(iout);
  jp = INTEGER_POINTER(jout);
  kp = INTEGER_POINTER(kout);
  xp = NUMERIC_POINTER(xout);
  yp = NUMERIC_POINTER(yout);
  sp = NUMERIC_POINTER(sinout);
  for(m = 0; m < Ne; m++) {
    ip[m] = ie[m];
    jp[m] = je[m];
    kp[m] = ke[m];
    xp[m] = xe[m];
    yp[m] = ye[m];
    sp[m] = sinalpha[m];
  }
  SET_VECTOR_ELT(out, 0, xout);
  SET_VECTOR_ELT(out, 1, yout);
  SET_VECTOR_ELT(out, 2, iout);
  SET_VECTOR_ELT(out, 3, jout);
  SET_VECTOR_ELT(out, 4, kout);
  SET_VECTOR_ELT(out, 5, sinout);    
  UNPROTECT(14);
  return(out);
}

/* circMseg: matrix of radii with rows corresponding to centres */

SEXP circMseg(SEXP XC,        /* circle centres */
              SEXP YC,        
	      SEXP R,         /* radii */
	      SEXP X0,        /* segments */
              SEXP Y0, 
              SEXP X1, 
	      SEXP Y1         
	      )
{
  double *xc, *yc, *r, *x0, *y0, *x1, *y1;
  int Nc, Ns, Nr, Ne, NeMax, newmax;

  /* outputs */
  int *ie, *je, *ke;  /* provenance of each intersection */
  double *xe, *ye, *sinalpha; /* cut points and angles */
  SEXP out, iout, jout, kout, xout, yout, sinout;
  int *ip, *jp, *kp;
  double *xp, *yp, *sp;

  /* internal */
  int i, j, k, m;
  double xci, yci, rik, x0c, y0c, dx, dy, A, B, C, Det, sqrtDet, sina;
  double u, u1, u2, slope, intcept, xcut, ycut, xnorm, ynorm, hx, hy;
  double twoA, fourA, Bsquared, Cdist2;

  PROTECT(XC = AS_NUMERIC(XC));
  PROTECT(YC = AS_NUMERIC(YC));
  PROTECT(R  = AS_NUMERIC(R));
  PROTECT(X0 = AS_NUMERIC(X0));
  PROTECT(Y0 = AS_NUMERIC(Y0));
  PROTECT(X1 = AS_NUMERIC(X1));
  PROTECT(Y1 = AS_NUMERIC(Y1));
  /* That's 7 protected */

  /* get pointers */
  xc = NUMERIC_POINTER(XC);
  yc = NUMERIC_POINTER(YC);
  r  = NUMERIC_POINTER(R);
  x0 = NUMERIC_POINTER(X0);
  y0 = NUMERIC_POINTER(Y0);
  x1 = NUMERIC_POINTER(X1);
  y1 = NUMERIC_POINTER(Y1);

  /* determine lengths of vectors */
  Nc = LENGTH(XC);
  Nr = LENGTH(R)/Nc;   /* n.b. */
  Ns = LENGTH(X0);

  /* Guess amount of storage required */
  NeMax = 4 * Nr * Nc;
  if(NeMax < 2048) NeMax = 2048;
  ie = (int *) R_alloc(NeMax, sizeof(int));
  je = (int *) R_alloc(NeMax, sizeof(int));
  ke = (int *) R_alloc(NeMax, sizeof(int));
  xe = (double *) R_alloc(NeMax, sizeof(double));
  ye = (double *) R_alloc(NeMax, sizeof(double));
  sinalpha = (double *) R_alloc(NeMax, sizeof(double));

  /* initialise output */
  Ne = 0;

  if(Nc > 0 && Ns > 0 && Nr > 0) {
    /* loop over circle centres */
    for(i = 0; i < Nc; i++) {
#ifdef BUGGER
      Rprintf("Circle %d\n", i);
#endif
      R_CheckUserInterrupt();
      xci = xc[i];
      yci = yc[i];
      /* loop over segments */
      for(j = 0; j < Ns; j++) {
#ifdef BUGGER
	Rprintf("\tSegment %d\n", j);
#endif
	dx = x1[j] - x0[j];
	dy = y1[j] - y0[j];
	x0c = x0[j] - xci;
	y0c = y0[j] - yci;
	/* find intersections between circle and infinite line */
	A = dx * dx + dy * dy;
	B = 2 * (dx * x0c + dy * y0c);
	twoA = 2.0 * A;
	fourA = 4.0 * A;
	Bsquared = B * B;
	Cdist2 = x0c * x0c + y0c * y0c;
	/* loop over radii */
	for(k = 0; k < Nr; k++) {
#ifdef BUGGER
	  Rprintf("\t\tRadius [%d, %d]\n", i, k);
#endif
	  rik = r[i + k*Nc];
	  C = Cdist2 - rik * rik;
	  Det = Bsquared - fourA * C;
	  if(Det > 0.0) {
	    /* two intersection points */
	    sqrtDet = sqrt(Det);
	    u1 = (-B - sqrtDet)/twoA;
	    u2 = (-B + sqrtDet)/twoA;
	    if(u1 > 0.0 && u1 < 1.0) {
	      /* first intersection point is inside segment */
	      if(dx != 0) {
		/* sloping line */
		slope = dy/dx;
		intcept = y0c - slope * x0c;
		xcut = x0c + u1 * dx;
		ycut = y0c + u1 * dy;
		ynorm = intcept/(slope * slope + 1.0);
		xnorm = - ynorm * slope;
	      } else {
		/* vertical line */
		xcut = x0c;
		ycut = y0c + u1 * dy;
		xnorm = x0c;
		ynorm = 0.0;
	      }
	      hx = xcut - xnorm;
	      hy = ycut - ynorm;

	      sina = sqrt(hx * hx + hy * hy)/rik;
	      if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	      /* add to output */
#ifdef BUGGER
	      Rprintf("\t\t\tAdding..\n");
#endif
	      sinalpha[Ne] = sina;
	      xe[Ne] = xcut + xci;
	      ye[Ne] = ycut + yci;
	      ie[Ne] = i + 1;
	      je[Ne] = j + 1;
	      ke[Ne] = k + 1;
	      ++Ne;
	      if(Ne >= NeMax) {
		/* storage overflow; reallocate */
#ifdef BUGGER
		Rprintf("\t\t\tOVERFLOW..\n");
#endif
		newmax = 2 * NeMax;
		xe = (double *) S_realloc((char *) xe,
					  newmax, NeMax, sizeof(double));
		ye = (double *) S_realloc((char *) ye,
					  newmax, NeMax, sizeof(double));
		ie = (int *) S_realloc((char *) ie,
				       newmax, NeMax, sizeof(int));
		je = (int *) S_realloc((char *) je,
				       newmax, NeMax, sizeof(int));
		ke = (int *) S_realloc((char *) ke,
				       newmax, NeMax, sizeof(int));
		sinalpha = (double *) S_realloc((char *) sinalpha,
						newmax, NeMax, sizeof(double));
		NeMax = newmax;
	      }
	    }
	    if(u2 > 0.0 && u2 < 1.0) {
	      /* second intersection point is inside segment */
	      if(dx != 0) {
		/* sloping line */
		slope = dy/dx;
		intcept = y0c - slope * x0c;
		xcut = x0c + u2 * dx;
		ycut = y0c + u2 * dy;
		ynorm = intcept/(slope * slope + 1.0);
		xnorm = - ynorm * slope;
	      } else {
		/* vertical line */
		xcut = x0c;
		ycut = y0c + u2 * dy;
		xnorm = x0c;
		ynorm = 0.0;
	      }
	      hx = xcut - xnorm;
	      hy = ycut - ynorm;

	      sina = sqrt(hx * hx + hy * hy)/rik;
	      if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	      /* add to output */
#ifdef BUGGER
	      Rprintf("\t\t\tAdding..\n");
#endif
	      sinalpha[Ne] = sina;
	      xe[Ne] = xcut + xci;
	      ye[Ne] = ycut + yci;
	      ie[Ne] = i + 1;
	      je[Ne] = j + 1;
	      ke[Ne] = k + 1;
	      ++Ne;
	      if(Ne >= NeMax) {
		/* storage overflow; reallocate */
#ifdef BUGGER
		Rprintf("\t\t\tOVERFLOW..\n");
#endif
		newmax = 2 * NeMax;
		xe = (double *) S_realloc((char *) xe,
					  newmax, NeMax, sizeof(double));
		ye = (double *) S_realloc((char *) ye,
					  newmax, NeMax, sizeof(double));
		ie = (int *) S_realloc((char *) ie,
				       newmax, NeMax, sizeof(int));
		je = (int *) S_realloc((char *) je,
				       newmax, NeMax, sizeof(int));
		ke = (int *) S_realloc((char *) ke,
				       newmax, NeMax, sizeof(int));
		sinalpha = (double *) S_realloc((char *) sinalpha,
						newmax, NeMax, sizeof(double));
		NeMax = newmax;
	      }
	    }
	  } else if(Det == 0.0) {
	    /* tangent point */
	    u = -B/twoA;
	    if(u > 0.0 && u < 1.0) {
	      /* tangent point is inside segment */
	      if(dx != 0) {
		/* sloping line */
		slope = dy/dx;
		intcept = y0c - slope * x0c;
		xcut = x0c + u * dx;
		ycut = y0c + u * dy;
		ynorm = intcept/(slope * slope + 1.0);
		xnorm = - ynorm * slope;
	      } else {
		/* vertical line */
		xcut = x0c;
		ycut = y0c + u * dy;
		xnorm = x0c;
		ynorm = 0.0;
	      }
	      hx = xcut - xnorm;
	      hy = ycut - ynorm;

	      sina = sqrt(hx * hx + hy * hy)/rik;
	      if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	      /* add to output */
#ifdef BUGGER
	      Rprintf("\t\t\tAdding..\n");
#endif
	      sinalpha[Ne] = sina;
	      xe[Ne] = xcut + xci;
	      ye[Ne] = ycut + yci;
	      ie[Ne] = i + 1;
	      je[Ne] = j + 1;
	      ke[Ne] = k + 1;
	      ++Ne;
	      if(Ne >= NeMax) {
		/* storage overflow; reallocate */
#ifdef BUGGER
		Rprintf("\t\t\tOVERFLOW..\n");
#endif
		newmax = 2 * NeMax;
		xe = (double *) S_realloc((char *) xe,
					  newmax, NeMax, sizeof(double));
		ye = (double *) S_realloc((char *) ye,
					  newmax, NeMax, sizeof(double));
		ie = (int *) S_realloc((char *) ie,
				       newmax, NeMax, sizeof(int));
		je = (int *) S_realloc((char *) je,
				       newmax, NeMax, sizeof(int));
		ke = (int *) S_realloc((char *) ke,
				       newmax, NeMax, sizeof(int));
		sinalpha = (double *) S_realloc((char *) sinalpha,
						newmax, NeMax, sizeof(double));
		NeMax = newmax;
	      }
	    }
	  }
	}
      }
    }
  }

  /* pack up */
  PROTECT(out = NEW_LIST(6));
  PROTECT(iout = NEW_INTEGER(Ne));
  PROTECT(jout = NEW_INTEGER(Ne));
  PROTECT(kout = NEW_INTEGER(Ne));
  PROTECT(xout = NEW_NUMERIC(Ne));  
  PROTECT(yout = NEW_NUMERIC(Ne));  
  PROTECT(sinout = NEW_NUMERIC(Ne));  
  /* 7 + 1 + 6 = 14 protected */
  ip = INTEGER_POINTER(iout);
  jp = INTEGER_POINTER(jout);
  kp = INTEGER_POINTER(kout);
  xp = NUMERIC_POINTER(xout);
  yp = NUMERIC_POINTER(yout);
  sp = NUMERIC_POINTER(sinout);
  for(m = 0; m < Ne; m++) {
    ip[m] = ie[m];
    jp[m] = je[m];
    kp[m] = ke[m];
    xp[m] = xe[m];
    yp[m] = ye[m];
    sp[m] = sinalpha[m];
  }
  SET_VECTOR_ELT(out, 0, xout);
  SET_VECTOR_ELT(out, 1, yout);
  SET_VECTOR_ELT(out, 2, iout);
  SET_VECTOR_ELT(out, 3, jout);
  SET_VECTOR_ELT(out, 4, kout);
  SET_VECTOR_ELT(out, 5, sinout);    
  UNPROTECT(14);
  return(out);
}

/* circPseg: centres and radii matched ('vectors') */

SEXP circPseg(SEXP XC,        /* circles (x0, y0, r) */
              SEXP YC,        
	      SEXP RC,         
	      SEXP X0,        /* segments */
              SEXP Y0, 
              SEXP X1, 
	      SEXP Y1         
	      )
{
  double *xc, *yc, *rc, *x0, *y0, *x1, *y1;
  int Nc, Ns, Ne, NeMax, newmax;

  /* outputs */
  int *ie, *je;  /* provenance of each intersection */
  double *xe, *ye, *sinalpha; /* cut points and angles */
  SEXP out, iout, jout, xout, yout, sinout;
  int *ip, *jp;
  double *xp, *yp, *sp;

  /* internal */
  int i, j, m;
  double xci, yci, rci, x0c, y0c, dx, dy, A, B, C, Det, sqrtDet, sina;
  double u, u1, u2, slope, intcept, xcut, ycut, xnorm, ynorm, hx, hy;
  double twoA, rci2;

  PROTECT(XC = AS_NUMERIC(XC));
  PROTECT(YC = AS_NUMERIC(YC));
  PROTECT(RC  = AS_NUMERIC(RC));
  PROTECT(X0 = AS_NUMERIC(X0));
  PROTECT(Y0 = AS_NUMERIC(Y0));
  PROTECT(X1 = AS_NUMERIC(X1));
  PROTECT(Y1 = AS_NUMERIC(Y1));
  /* That's 7 protected */

  /* get pointers */
  xc = NUMERIC_POINTER(XC);
  yc = NUMERIC_POINTER(YC);
  rc = NUMERIC_POINTER(RC);
  x0 = NUMERIC_POINTER(X0);
  y0 = NUMERIC_POINTER(Y0);
  x1 = NUMERIC_POINTER(X1);
  y1 = NUMERIC_POINTER(Y1);

  /* determine lengths of vectors */
  Nc = LENGTH(XC);
  Ns = LENGTH(X0);

  /* Guess amount of storage required */
  NeMax = 4 * (Ns + Nc);
  if(NeMax < 2048) NeMax = 2048;
  ie = (int *) R_alloc(NeMax, sizeof(int));
  je = (int *) R_alloc(NeMax, sizeof(int));
  xe = (double *) R_alloc(NeMax, sizeof(double));
  ye = (double *) R_alloc(NeMax, sizeof(double));
  sinalpha = (double *) R_alloc(NeMax, sizeof(double));

  /* initialise output */
  Ne = 0;

  if(Nc > 0 && Ns > 0) {
    /* loop over circles */
    for(i = 0; i < Nc; i++) {

#ifdef BUGGER
      Rprintf("Circle %d\n", i);
#endif
      R_CheckUserInterrupt();

      xci = xc[i];
      yci = yc[i];
      rci = rc[i];

      rci2 = rci * rci;

      /* loop over segments */
      for(j = 0; j < Ns; j++) {

#ifdef BUGGER
	Rprintf("\tSegment %d\n", j);
#endif

	dx = x1[j] - x0[j];
	dy = y1[j] - y0[j];
	x0c = x0[j] - xci;
	y0c = y0[j] - yci;
	/* find intersections between circle and infinite line */
	A = dx * dx + dy * dy;
	B = 2 * (dx * x0c + dy * y0c);
	C = x0c * x0c + y0c * y0c - rci2;
	Det = B * B - 4.0 * A * C;
	twoA = 2.0 * A;
	if(Det > 0.0) {
	  /* two intersection points */
	  sqrtDet = sqrt(Det);
	  u1 = (-B - sqrtDet)/twoA;
	  u2 = (-B + sqrtDet)/twoA;
	  if(u1 > 0.0 && u1 < 1.0) {
	    /* first intersection point is inside segment */
	    if(dx != 0) {
	      /* sloping line */
	      slope = dy/dx;
	      intcept = y0c - slope * x0c;
	      xcut = x0c + u1 * dx;
	      ycut = y0c + u1 * dy;
	      ynorm = intcept/(slope * slope + 1.0);
	      xnorm = - ynorm * slope;
	    } else {
	      /* vertical line */
	      xcut = x0c;
	      ycut = y0c + u1 * dy;
	      xnorm = x0c;
	      ynorm = 0.0;
	    }
	    hx = xcut - xnorm;
	    hy = ycut - ynorm;
	    
	    sina = sqrt(hx * hx + hy * hy)/rci;
	    if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	    /* add to output */
#ifdef BUGGER
	    Rprintf("\t\t\tAdding..\n");
#endif
	    sinalpha[Ne] = sina;
	    xe[Ne] = xcut + xci;
	    ye[Ne] = ycut + yci;
	    ie[Ne] = i + 1;
	    je[Ne] = j + 1;
	    ++Ne;
	    if(Ne >= NeMax) {
	      /* storage overflow; reallocate */
#ifdef BUGGER
	      Rprintf("\t\t\tOVERFLOW..\n");
#endif
	      newmax = 2 * NeMax;
	      xe = (double *) S_realloc((char *) xe,
					newmax, NeMax, sizeof(double));
	      ye = (double *) S_realloc((char *) ye,
					newmax, NeMax, sizeof(double));
	      ie = (int *) S_realloc((char *) ie,
				     newmax, NeMax, sizeof(int));
	      je = (int *) S_realloc((char *) je,
				     newmax, NeMax, sizeof(int));
	      sinalpha = (double *) S_realloc((char *) sinalpha,
					      newmax, NeMax, sizeof(double));
	      NeMax = newmax;
	    }
	  }
	  if(u2 > 0.0 && u2 < 1.0) {
	    /* second intersection point is inside segment */
	    if(dx != 0) {
	      /* sloping line */
	      slope = dy/dx;
	      intcept = y0c - slope * x0c;
	      xcut = x0c + u2 * dx;
	      ycut = y0c + u2 * dy;
	      ynorm = intcept/(slope * slope + 1.0);
	      xnorm = - ynorm * slope;
	    } else {
	      /* vertical line */
	      xcut = x0c;
	      ycut = y0c + u2 * dy;
	      xnorm = x0c;
	      ynorm = 0.0;
	    }
	    hx = xcut - xnorm;
	    hy = ycut - ynorm;
	    
	    sina = sqrt(hx * hx + hy * hy)/rci;
	    if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	    /* add to output */
#ifdef BUGGER
	    Rprintf("\t\t\tAdding..\n");
#endif
	    sinalpha[Ne] = sina;
	    xe[Ne] = xcut + xci;
	    ye[Ne] = ycut + yci;
	    ie[Ne] = i + 1;
	    je[Ne] = j + 1;
	    ++Ne;
	    if(Ne >= NeMax) {
	      /* storage overflow; reallocate */
#ifdef BUGGER
	      Rprintf("\t\t\tOVERFLOW..\n");
#endif
	      newmax = 2 * NeMax;
	      xe = (double *) S_realloc((char *) xe,
					newmax, NeMax, sizeof(double));
	      ye = (double *) S_realloc((char *) ye,
					newmax, NeMax, sizeof(double));
	      ie = (int *) S_realloc((char *) ie,
				     newmax, NeMax, sizeof(int));
	      je = (int *) S_realloc((char *) je,
				     newmax, NeMax, sizeof(int));
	      sinalpha = (double *) S_realloc((char *) sinalpha,
					      newmax, NeMax, sizeof(double));
	      NeMax = newmax;
	    }
	  }
	} else if(Det == 0.0) {
	  /* tangent point */
	  u = -B/twoA;
	  if(u > 0.0 && u < 1.0) {
	    /* tangent point is inside segment */
	    if(dx != 0) {
	      /* sloping line */
	      slope = dy/dx;
	      intcept = y0c - slope * x0c;
	      xcut = x0c + u * dx;
	      ycut = y0c + u * dy;
	      ynorm = intcept/(slope * slope + 1.0);
	      xnorm = - ynorm * slope;
	    } else {
	      /* vertical line */
	      xcut = x0c;
	      ycut = y0c + u * dy;
	      xnorm = x0c;
	      ynorm = 0.0;
	    }
	    hx = xcut - xnorm;
	    hy = ycut - ynorm;
	    
	    sina = sqrt(hx * hx + hy * hy)/rci;
	    if(sina > 1.0) sina = 1.0; else if(sina < -1.0) sina = -1.0;

	    /* add to output */
#ifdef BUGGER
	    Rprintf("\t\t\tAdding..\n");
#endif
	    sinalpha[Ne] = sina;
	    xe[Ne] = xcut + xci;
	    ye[Ne] = ycut + yci;
	    ie[Ne] = i + 1;
	    je[Ne] = j + 1;
	    ++Ne;
	    if(Ne >= NeMax) {
	      /* storage overflow; reallocate */
#ifdef BUGGER
	      Rprintf("\t\t\tOVERFLOW..\n");
#endif
	      newmax = 2 * NeMax;
	      xe = (double *) S_realloc((char *) xe,
					newmax, NeMax, sizeof(double));
	      ye = (double *) S_realloc((char *) ye,
					newmax, NeMax, sizeof(double));
	      ie = (int *) S_realloc((char *) ie,
				     newmax, NeMax, sizeof(int));
	      je = (int *) S_realloc((char *) je,
				     newmax, NeMax, sizeof(int));
	      sinalpha = (double *) S_realloc((char *) sinalpha,
					      newmax, NeMax, sizeof(double));
	      NeMax = newmax;
	    }
	  }
	}
      }
    }
  }

  /* pack up */
  PROTECT(out = NEW_LIST(5));
  PROTECT(iout = NEW_INTEGER(Ne));
  PROTECT(jout = NEW_INTEGER(Ne));
  PROTECT(xout = NEW_NUMERIC(Ne));  
  PROTECT(yout = NEW_NUMERIC(Ne));  
  PROTECT(sinout = NEW_NUMERIC(Ne));  
  /* 7 + 1 + 5 = 13 protected */
  ip = INTEGER_POINTER(iout);
  jp = INTEGER_POINTER(jout);
  xp = NUMERIC_POINTER(xout);
  yp = NUMERIC_POINTER(yout);
  sp = NUMERIC_POINTER(sinout);
  for(m = 0; m < Ne; m++) {
    ip[m] = ie[m];
    jp[m] = je[m];
    xp[m] = xe[m];
    yp[m] = ye[m];
    sp[m] = sinalpha[m];
  }
  SET_VECTOR_ELT(out, 0, xout);
  SET_VECTOR_ELT(out, 1, yout);
  SET_VECTOR_ELT(out, 2, iout);
  SET_VECTOR_ELT(out, 3, jout);
  SET_VECTOR_ELT(out, 4, sinout);    
  UNPROTECT(13);
  return(out);
}
