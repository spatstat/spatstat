/*

  xyseg.c

  Computation with line segments

  xysegint     compute intersections between line segments

  $Revision: 1.20 $     $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "chunkloop.h"

#define NIETS -1.0

#undef DEBUG 
#define INSIDE01(X,E) (X * (1.0 - X) >= -E)

/* 
         --------------- PAIRS OF PSP OBJECTS ----------------------
*/

/*  
   xysegint

   Determines intersections between each pair of line segments
   drawn from two lists of line segments.

   Line segments are given as x0, y0, dx, dy
   where (x0,y0) is the first endpoint and (dx, dy) is the vector
   from the first to the second endpoint.
   Points along a line segment are represented in parametric
   coordinates, 
            (x,y) = (x0, y0) + t * (dx, dy).

   Output from xysegint() consists of five matrices xx, yy, ta, tb, ok.
   The (i,j)-th entries in these matrices give information about the
   intersection between the i-th segment in list 'a' and the
   j-th segment in list 'b'. The information is

       ok[i,j]  = 1 if there is an intersection
                = 0 if not

       xx[i,j]  = x coordinate of intersection

       yy[i,j]  = y coordinate of intersection

       ta[i,j] = parameter of intersection point
                 relative to i-th segment in list 'a'

       tb[i,j] = parameter of intersection point
                 relative to j-th segment in list 'b'

*/
	     

void xysegint(na, x0a, y0a, dxa, dya, 
              nb, x0b, y0b, dxb, dyb, 
	      eps,
              xx, yy, ta, tb, ok)
     /* inputs (vectors of coordinates) */
     int *na, *nb;
     double *x0a, *y0a, *dxa, *dya, *x0b, *y0b, *dxb, *dyb;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     double *xx, *yy, *ta, *tb;
     int *ok;
{ 
  int i, j, ma, mb, ijpos, maxchunk;
  double determinant, absdet, diffx, diffy, tta, ttb, epsilon;

  ma = *na;
  mb = *nb;
  epsilon = *eps;

  OUTERCHUNKLOOP(j, mb, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, mb, maxchunk, 8196) {
      for(i = 0; i < ma; i++) {
	ijpos = j * ma + i;
	ok[ijpos] = 0;
	xx[ijpos] = yy[ijpos] = ta[ijpos] = tb[ijpos] = NIETS;
	determinant = dxb[j] * dya[i] - dyb[j] * dxa[i];
	absdet = (determinant > 0) ? determinant : -determinant;
#ifdef DEBUG
	Rprintf("i = %d, j = %d\n", i, j);
	Rprintf("segment A[i]: (%lf, %lf) to (%lf, %lf)\n",
		x0a[i], y0a[i], x0a[i] + dxa[i], y0a[i] + dya[i]);
	Rprintf("segment B[j]: (%lf, %lf) to (%lf, %lf)\n",
		x0b[j], y0b[j], x0b[j] + dxb[j], y0b[j] + dyb[j]);
	Rprintf("determinant=%lf\n", determinant);
#endif	
	if(absdet > epsilon) {
	  diffx = (x0b[j] - x0a[i])/determinant;
	  diffy = (y0b[j] - y0a[i])/determinant;
	  ta[ijpos] = tta = - dyb[j] * diffx + dxb[j] * diffy;
	  tb[ijpos] = ttb = - dya[i] * diffx + dxa[i] * diffy;
#ifdef DEBUG
	  Rprintf("ta = %lf, tb = %lf\n", tta, ttb);
#endif	
	  if(INSIDE01(tta, epsilon) && INSIDE01(ttb, epsilon)) {
	    /* intersection */
	    ok[ijpos] = 1;
	    xx[ijpos] = x0a[i] + tta * dxa[i];
	    yy[ijpos] = y0a[i] + tta * dya[i];
#ifdef DEBUG
	    Rprintf("segments intersect at (%lf, %lf)\n", xx[ijpos], yy[ijpos]);
#endif	
	  }
	}
      }
    }
  }
}

/* 
   Stripped-down version of xysegint that just returns logical matrix 
*/

void xysi(na, x0a, y0a, dxa, dya, 
              nb, x0b, y0b, dxb, dyb, 
	      eps,
              ok)
     /* inputs (vectors of coordinates) */
     int *na, *nb;
     double *x0a, *y0a, *dxa, *dya, *x0b, *y0b, *dxb, *dyb;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     int *ok;
{ 
  int i, j, ma, mb, ijpos, maxchunk;
  double determinant, absdet, diffx, diffy, tta, ttb, epsilon;

  ma = *na;
  mb = *nb;
  epsilon = *eps;

  OUTERCHUNKLOOP(j, mb, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, mb, maxchunk, 8196) {
      for(i = 0; i < ma; i++) {
	ijpos = j * ma + i;
	ok[ijpos] = 0;
	determinant = dxb[j] * dya[i] - dyb[j] * dxa[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > epsilon) {
	  diffx = (x0b[j] - x0a[i])/determinant;
	  diffy = (y0b[j] - y0a[i])/determinant;
	  tta = - dyb[j] * diffx + dxb[j] * diffy;
	  ttb = - dya[i] * diffx + dxa[i] * diffy;
	  if(INSIDE01(tta, epsilon) && INSIDE01(ttb, epsilon)) {
	    /* intersection */
	    ok[ijpos] = 1;
	  }
	}
      }
    }
  }
}

/* 
   Test whether there is at least one intersection
*/

void xysiANY(na, x0a, y0a, dxa, dya, 
		nb, x0b, y0b, dxb, dyb, 
		eps,
		ok)
     /* inputs (vectors of coordinates) */
     int *na, *nb;
     double *x0a, *y0a, *dxa, *dya, *x0b, *y0b, *dxb, *dyb;
     /* input (tolerance for determinant) */
     double *eps;  
     /* output (single logical value) */
     int *ok;
{ 
  int i, j, ma, mb, maxchunk;
  double determinant, absdet, diffx, diffy, tta, ttb, epsilon;

  *ok = 0;
  ma = *na;
  mb = *nb;
  epsilon = *eps;

  OUTERCHUNKLOOP(j, mb, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, mb, maxchunk, 8196) {
      for(i = 0; i < ma; i++) {
	determinant = dxb[j] * dya[i] - dyb[j] * dxa[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > epsilon) {
	  diffx = (x0b[j] - x0a[i])/determinant;
	  diffy = (y0b[j] - y0a[i])/determinant;
	  tta = - dyb[j] * diffx + dxb[j] * diffy;
	  ttb = - dya[i] * diffx + dxa[i] * diffy;
	  if(INSIDE01(tta, epsilon) && INSIDE01(ttb, epsilon)) {
	    /* intersection */
	    *ok = 1;
	    return;
	  }
	}
      }
    }
  }
}

/* 
    Analogue of xysegint
    when segments in list 'a' are infinite vertical lines
*/

void xysegVslice(na, xa,  
		 nb, x0b, y0b, dxb, dyb, 
		 eps,
		 yy, ok)
     /* inputs (vectors of coordinates) */
     int *na, *nb;
     double *xa, *x0b, *y0b, *dxb, *dyb;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     double *yy;
     int *ok;
{ 
  int i, j, ma, mb, ijpos, maxchunk;
  double diffx0, diffx1, width, abswidth, epsilon;
  int notvertical;

  ma = *na;
  mb = *nb;
  epsilon = *eps;

  OUTERCHUNKLOOP(j, mb, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, mb, maxchunk, 8196) {
      /* determine whether segment j is nearly vertical */
      width = dxb[j];
      abswidth = (width > 0) ? width : -width;
      notvertical = (abswidth <= epsilon);
    
      for(i = 0; i < ma; i++) {
	ijpos = j * ma + i;
	ok[ijpos] = 0;
	yy[ijpos] = NIETS;
	/* test whether vertical line i separates endpoints of segment j */
	diffx0 = xa[i] - x0b[j];
	diffx1 = diffx0 - width;
	if(diffx0 * diffx1 <= 0) {
	  /* intersection */
	  ok[ijpos] = 1;
	  /* compute y-coordinate of intersection point */
	  if(notvertical) {
	    yy[ijpos] = y0b[j] + diffx0 * dyb[j]/width;
	  } else {
	    /* vertical or nearly-vertical segment: pick midpoint */	  
	    yy[ijpos] = y0b[j] + dyb[j]/2.0;
	  }
	}
      }
    }
  }
}

/* 
    -------------- ONE PSP OBJECT ----------------------------
*/
	 

/*

    Similar to xysegint,
    but computes intersections between all pairs of segments
    in a single list, excluding the diagonal comparisons of course

*/

void xysegXint(n, x0, y0, dx, dy, 
	      eps,
              xx, yy, ti, tj, ok)
     /* inputs (vectors of coordinates) */
     int *n;
     double *x0, *y0, *dx, *dy;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     double *xx, *yy, *ti, *tj;
     int *ok;
{ 
  int i, j, m, mm1, ijpos, jipos, iipos, maxchunk;
  double determinant, absdet, diffx, diffy, tti, ttj, epsilon;

  m = *n;
  epsilon = *eps;
 
  mm1 = m - 1;
  OUTERCHUNKLOOP(j, mm1, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, mm1, maxchunk, 8196) {
      for(i = j+1; i < m; i++) {
	ijpos = j * m + i;
	jipos = i * m + j;
	ok[ijpos] = ok[jipos] = 0;
	xx[ijpos] = yy[ijpos] = ti[ijpos] = ti[jipos] = NIETS;
	xx[jipos] = yy[jipos] = tj[ijpos] = tj[jipos] = NIETS;
	determinant = dx[j] * dy[i] - dy[j] * dx[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > epsilon) {
	  diffx = (x0[j] - x0[i])/determinant;
	  diffy = (y0[j] - y0[i])/determinant;
	  ti[ijpos] = tti = - dy[j] * diffx + dx[j] * diffy;
	  tj[ijpos] = ttj = - dy[i] * diffx + dx[i] * diffy;
	  tj[jipos] = ti[ijpos];
	  ti[jipos] = tj[ijpos];
	  if(INSIDE01(tti, epsilon) && INSIDE01(ttj, epsilon)) {
	    ok[ijpos] = ok[jipos] = 1;
	    xx[ijpos] = xx[jipos] = x0[i] + tti * dx[i];
	    yy[ijpos] = yy[jipos] = y0[i] + tti * dy[i];
	  }
	}
      }
    }
  }

  /* assign diagonal */
  for(i = 0; i < m; i++) {
    iipos = i * m + i;
    ok[iipos] = 0;
    xx[iipos] = yy[iipos] = ti[iipos] = tj[iipos] = NIETS;
  }

}
	 
/*

    Reduced version of xysegXint that returns logical matrix 'ok' only

*/

void xysxi(n, x0, y0, dx, dy, 
	      eps,
              ok)
     /* inputs (vectors of coordinates) */
     int *n;
     double *x0, *y0, *dx, *dy;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     int *ok;
{ 
  int i, j, m, mm1, ijpos, jipos, iipos, maxchunk;
  double determinant, absdet, diffx, diffy, tti, ttj, epsilon;

  m = *n;
  epsilon = *eps;
 
  mm1 = m - 1;
  OUTERCHUNKLOOP(j, mm1, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, mm1, maxchunk, 8196) {
      for(i = j+1; i < m; i++) {
	ijpos = j * m + i;
	jipos = i * m + j;
	ok[ijpos] = ok[jipos] = 0;
	determinant = dx[j] * dy[i] - dy[j] * dx[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > epsilon) {
	  diffx = (x0[j] - x0[i])/determinant;
	  diffy = (y0[j] - y0[i])/determinant;
	  tti = - dy[j] * diffx + dx[j] * diffy;
	  ttj = - dy[i] * diffx + dx[i] * diffy;
	  if(INSIDE01(tti, epsilon) && INSIDE01(ttj, epsilon)) {
	    ok[ijpos] = ok[jipos] = 1;
	  }
	}
      }
    }
  }

  /* assign diagonal */
  for(i = 0; i < m; i++) {
    iipos = i * m + i;
    ok[iipos] = 0;
  }

}

/*
   ---------------------- ONE CLOSED POLYGON ------------------------
*/
	 
/*

    Identify self-intersections in a closed polygon

    (Similar to xysegXint,
    but does not compare segments which are cyclically adjacent in the list)

*/

void Cxypolyselfint(n, x0, y0, dx, dy, 
	      eps,
              xx, yy, ti, tj, ok)
     /* inputs (vectors of coordinates) */
     int *n;
     double *x0, *y0, *dx, *dy;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     double *xx, *yy, *ti, *tj;
     int *ok;
{ 
  int i, j, k, m, m2, mm1, mm2, mstop, ijpos, jipos, maxchunk;
  double determinant, absdet, diffx, diffy, tti, ttj, epsilon;

  m = *n;
  epsilon = *eps;
  m2 = m * m;

  /* initialise matrices */
  
  for(k = 0; k < m2; k++) {
    ok[k] = 0;
    xx[k] = yy[k] = ti[k] = tj[k] = NIETS;
  }

  if(m <= 2) 
    return;

  /* Compare j with j+2, j+3, ...., m-1
     Don't compare 0 with m-1  */
  mm1 = m - 1;
  mm2 = m - 2;
  OUTERCHUNKLOOP(j, mm2, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, mm2, maxchunk, 8196) {
      mstop = (j > 0) ? m : mm1;
      for(i = j+2; i < mstop; i++) {
	ijpos = j * m + i;
	jipos = i * m + j;
	determinant = dx[j] * dy[i] - dy[j] * dx[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > epsilon) {
	  diffx = (x0[j] - x0[i])/determinant;
	  diffy = (y0[j] - y0[i])/determinant;
	  ti[ijpos] = tti = - dy[j] * diffx + dx[j] * diffy;
	  tj[ijpos] = ttj = - dy[i] * diffx + dx[i] * diffy;
	  tj[jipos] = ti[ijpos];
	  ti[jipos] = tj[ijpos];
	  if(INSIDE01(tti, epsilon) && INSIDE01(ttj, epsilon)) {
	    ok[ijpos] = ok[jipos] = 1;
	    xx[ijpos] = xx[jipos] = x0[i] + tti * dx[i];
	    yy[ijpos] = yy[jipos] = y0[i] + tti * dy[i];
	  }
	}
      }
    }
  }
}
	 

/*
  Just determines whether there is self-intersection
  (exits quicker & uses less space)
*/


void xypsi(n, x0, y0, dx, dy, xsep, ysep, eps, proper, answer)
     /* inputs (vectors of coordinates) */
     int *n;
     double *x0, *y0, *dx, *dy;
     /* inputs (distances beyond which intersection is impossible) */
     double *xsep, *ysep;
     /* input (tolerance for determinant) */
     double *eps;  
     /* input (flag) */
     int *proper;
     /* output */
     int *answer;
{ 
  int i, j, m, mm1, mm2, mstop, prop, maxchunk;
  double determinant, absdet, diffx, diffy, tti, ttj, epsilon;
  double Xsep, Ysep;

  m = *n;
  prop = *proper;
  Xsep = *xsep;
  Ysep = *ysep;
  epsilon = *eps;

  *answer = 0;

  if(m <= 2) 
    return;

  /* Compare j with j+2, j+3, ...., m-1
     Don't compare 0 with m-1  */
  mm1 = m - 1;
  mm2 = m - 2;
  OUTERCHUNKLOOP(j, mm2, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, mm2, maxchunk, 8196) {
      mstop = (j > 0) ? m : mm1;
      for(i = j+2; i < mstop; i++) {
	diffx = x0[j] - x0[i];
	diffy = y0[j] - y0[i];
	if(diffx < Xsep && diffx > -Xsep && diffy < Ysep && diffy > -Ysep) {
	  determinant = dx[j] * dy[i] - dy[j] * dx[i];
	  absdet = (determinant > 0) ? determinant : -determinant;
	  if(absdet > epsilon) {
	    diffx = diffx/determinant;
	    diffy = diffy/determinant;
	    tti = - dy[j] * diffx + dx[j] * diffy;
	    ttj = - dy[i] * diffx + dx[i] * diffy;
	    if(INSIDE01(tti, epsilon) && INSIDE01(ttj, epsilon)) {
              /* intersection occurs */
	      if(prop == 0 ||
		 (tti != 0.0 && tti != 1.0) || 
		 (ttj != 0.0 && ttj != 1.0)) {
              /* proper intersection */
		*answer = 1;
		return;
	      }
	    }
	  }
	}
      }
    }
  }
}

	 
/*
        ---------------- .Call INTERFACE  ---------------------------

	Analogues of functions above, but using the .Call interface
	and dynamic storage allocation, to save space.

 */

SEXP Cxysegint(SEXP x0a, 
               SEXP y0a, 
               SEXP dxa, 
               SEXP dya, 
               SEXP x0b, 
               SEXP y0b, 
               SEXP dxb, 
               SEXP dyb, 
	       SEXP eps) 
{
  int i, j, k, na, nb;
  double determinant, absdet, diffx, diffy, tta, ttb;

  int nout, noutmax, newmax, maxchunk;
  double epsilon;
  double *x0A, *y0A, *dxA, *dyA, *x0B, *y0B, *dxB, *dyB;
  double *ta, *tb, *x, *y;
  int *ia, *jb;
  SEXP out, iAout, jBout, tAout, tBout, xout, yout;
  double *tAoutP, *tBoutP, *xoutP, *youtP;
  int *iAoutP, *jBoutP;

  PROTECT(x0a = AS_NUMERIC(x0a));
  PROTECT(y0a = AS_NUMERIC(y0a));
  PROTECT(dxa = AS_NUMERIC(dxa));
  PROTECT(dya = AS_NUMERIC(dya));
  PROTECT(x0b = AS_NUMERIC(x0b));
  PROTECT(y0b = AS_NUMERIC(y0b));
  PROTECT(dxb = AS_NUMERIC(dxb));
  PROTECT(dyb = AS_NUMERIC(dyb));
  PROTECT(eps = AS_NUMERIC(eps));
  /* that's 9 protected */


  /* get pointers */
  x0A = NUMERIC_POINTER(x0a);
  y0A = NUMERIC_POINTER(y0a);
  dxA = NUMERIC_POINTER(dxa);
  dyA = NUMERIC_POINTER(dya);
  x0B = NUMERIC_POINTER(x0b);
  y0B = NUMERIC_POINTER(y0b);
  dxB = NUMERIC_POINTER(dxb);
  dyB = NUMERIC_POINTER(dyb);

  /* determine length of vectors */
  na = LENGTH(x0a);
  nb = LENGTH(x0b);
  epsilon = *(NUMERIC_POINTER(eps));
  
  /* guess amount of storage required for output */
  noutmax = (na > nb) ? na : nb;
  nout = 0;
  ia = (int *) R_alloc(noutmax, sizeof(int));
  jb = (int *) R_alloc(noutmax, sizeof(int));
  ta = (double *) R_alloc(noutmax, sizeof(double));
  tb = (double *) R_alloc(noutmax, sizeof(double));
  x = (double *) R_alloc(noutmax, sizeof(double));
  y = (double *) R_alloc(noutmax, sizeof(double));

  /* scan data and collect intersections */
  OUTERCHUNKLOOP(j, nb, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nb, maxchunk, 8196) {
      for(i = 0; i < na; i++) {
	determinant = dxB[j] * dyA[i] - dyB[j] * dxA[i];
	absdet = (determinant > 0) ? determinant : -determinant;
#ifdef DEBUG
	Rprintf("i = %d, j = %d\n", i, j);
	Rprintf("segment A[i]: (%lf, %lf) to (%lf, %lf)\n",
		x0A[i], y0A[i], x0A[i] + dxA[i], y0A[i] + dyA[i]);
	Rprintf("segment B[j]: (%lf, %lf) to (%lf, %lf)\n",
		x0B[j], y0B[j], x0B[j] + dxB[j], y0B[j] + dyB[j]);
	Rprintf("determinant=%lf\n", determinant);
#endif	
	if(absdet > epsilon) {
	  diffx = (x0B[j] - x0A[i])/determinant;
	  diffy = (y0B[j] - y0A[i])/determinant;
	  tta = - dyB[j] * diffx + dxB[j] * diffy;
	  ttb = - dyA[i] * diffx + dxA[i] * diffy;
#ifdef DEBUG
	  Rprintf("ta = %lf, tb = %lf\n", tta, ttb);
#endif	
	  if(INSIDE01(tta, epsilon) && INSIDE01(ttb, epsilon)) {
	    /* intersection */
	    if(nout >= noutmax) {
	      /* storage overflow - increase space */
	      newmax = 4 * noutmax;
	      ia = (int *) S_realloc((char *) ia, 
				     newmax, noutmax, sizeof(int));
	      jb = (int *) S_realloc((char *) jb, 
				     newmax, noutmax, sizeof(int));
	      ta = (double *) S_realloc((char *) ta, 
					newmax, noutmax, sizeof(double));
	      tb = (double *) S_realloc((char *) tb, 
					newmax, noutmax, sizeof(double));
	      x = (double *) S_realloc((char *) x, 
				       newmax, noutmax, sizeof(double));
	      y = (double *) S_realloc((char *) y, 
				       newmax, noutmax, sizeof(double));
	      noutmax = newmax;
	    }
	    ta[nout] = tta;
	    tb[nout] = ttb;
	    ia[nout] = i;
	    jb[nout] = j;
	    x[nout]  = x0A[i] + tta * dxA[i];
	    y[nout]  = y0A[i] + tta * dyA[i];
#ifdef DEBUG
	    Rprintf("segments intersect at (%lf, %lf)\n", x[nout], y[nout]);
#endif	
	    ++nout;
	  }
	}
      }
    }
  }
  /* pack up */
  PROTECT(iAout = NEW_INTEGER(nout));
  PROTECT(jBout = NEW_INTEGER(nout));
  PROTECT(tAout = NEW_NUMERIC(nout));
  PROTECT(tBout = NEW_NUMERIC(nout));
  PROTECT(xout = NEW_NUMERIC(nout));
  PROTECT(yout = NEW_NUMERIC(nout));
  /* 9 + 6 = 15 protected */
  iAoutP = INTEGER_POINTER(iAout);
  jBoutP = INTEGER_POINTER(jBout);
  tAoutP = NUMERIC_POINTER(tAout);
  tBoutP = NUMERIC_POINTER(tBout);
  xoutP = NUMERIC_POINTER(xout);
  youtP = NUMERIC_POINTER(yout);
  for(k = 0; k < nout; k++) {
    iAoutP[k] = ia[k];
    jBoutP[k] = jb[k];
    tAoutP[k] = ta[k];
    tBoutP[k] = tb[k];
    xoutP[k]  = x[k];
    youtP[k]  = y[k];
  }
  PROTECT(out = NEW_LIST(6));
  /* 15 + 1 = 16 protected */
  SET_VECTOR_ELT(out, 0, iAout);
  SET_VECTOR_ELT(out, 1, jBout);
  SET_VECTOR_ELT(out, 2, tAout);
  SET_VECTOR_ELT(out, 3, tBout);
  SET_VECTOR_ELT(out, 4, xout);
  SET_VECTOR_ELT(out, 5, yout);
  UNPROTECT(16);
  return(out);
}


/*

    Similar to Cxysegint,
    but computes intersections between all pairs of segments
    in a single list, excluding the diagonal comparisons of course

*/

SEXP CxysegXint(SEXP x0, 
		SEXP y0, 
		SEXP dx, 
		SEXP dy,
		SEXP eps)
{ 
  int i, j, k, n, n1;
  double determinant, absdet, diffx, diffy, tti, ttj;

  int nout, noutmax, newmax, maxchunk;
  double epsilon;
  double *X0, *Y0, *Dx, *Dy;
  double *ti, *tj, *x, *y;
  int *ii, *jj;
  SEXP out, iout, jout, tiout, tjout, xout, yout;
  double *tioutP, *tjoutP, *xoutP, *youtP;
  int *ioutP, *joutP;

  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(y0 = AS_NUMERIC(y0));
  PROTECT(dx = AS_NUMERIC(dx));
  PROTECT(dy = AS_NUMERIC(dy));
  PROTECT(eps = AS_NUMERIC(eps));
  /* that's 5 protected */

  /* get pointers */
  X0 = NUMERIC_POINTER(x0);
  Y0 = NUMERIC_POINTER(y0);
  Dx = NUMERIC_POINTER(dx);
  Dy = NUMERIC_POINTER(dy);
  
  /* determine length of vectors */
  n = LENGTH(x0);
  epsilon = *(NUMERIC_POINTER(eps));

  /* guess amount of storage required for output */
  noutmax = n;
  nout = 0;
  ii = (int *) R_alloc(noutmax, sizeof(int));
  jj = (int *) R_alloc(noutmax, sizeof(int));
  ti = (double *) R_alloc(noutmax, sizeof(double));
  tj = (double *) R_alloc(noutmax, sizeof(double));
  x = (double *) R_alloc(noutmax, sizeof(double));
  y = (double *) R_alloc(noutmax, sizeof(double));

  /* scan data */
  n1 = n - 1;
  OUTERCHUNKLOOP(j, n1, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, n1, maxchunk, 8196) {
      for(i = j+1; i < n; i++) {
	determinant = Dx[j] * Dy[i] - Dy[j] * Dx[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > epsilon) {
	  diffx = (X0[j] - X0[i])/determinant;
	  diffy = (Y0[j] - Y0[i])/determinant;
	  tti = - Dy[j] * diffx + Dx[j] * diffy;
	  ttj = - Dy[i] * diffx + Dx[i] * diffy;
	  if(INSIDE01(tti,epsilon) && INSIDE01(ttj,epsilon)) {
	    /* intersection */
	    if(nout >= noutmax) {
	      /* storage overflow - increase space */
	      newmax = 4 * noutmax;
	      ii = (int *) S_realloc((char *) ii, 
				     newmax, noutmax, sizeof(int));
	      jj = (int *) S_realloc((char *) jj, 
				     newmax, noutmax, sizeof(int));
	      ti = (double *) S_realloc((char *) ti, 
					newmax, noutmax, sizeof(double));
	      tj = (double *) S_realloc((char *) tj, 
					newmax, noutmax, sizeof(double));
	      x = (double *) S_realloc((char *) x, 
				       newmax, noutmax, sizeof(double));
	      y = (double *) S_realloc((char *) y, 
				       newmax, noutmax, sizeof(double));
	      noutmax = newmax;
	    }
	    ti[nout] = tti;
	    tj[nout] = ttj;
	    ii[nout] = i;
	    jj[nout] = j;
	    x[nout]  = X0[i] + tti * Dx[i];
	    y[nout]  = Y0[i] + tti * Dy[i];
	    ++nout;
	  }
	}
      }
    }
  }

  /* pack up */
  PROTECT(iout = NEW_INTEGER(nout));
  PROTECT(jout = NEW_INTEGER(nout));
  PROTECT(tiout = NEW_NUMERIC(nout));
  PROTECT(tjout = NEW_NUMERIC(nout));
  PROTECT(xout = NEW_NUMERIC(nout));
  PROTECT(yout = NEW_NUMERIC(nout));
  /* 5 + 6 = 11 protected */
  ioutP = INTEGER_POINTER(iout);
  joutP = INTEGER_POINTER(jout);
  tioutP = NUMERIC_POINTER(tiout);
  tjoutP = NUMERIC_POINTER(tjout);
  xoutP = NUMERIC_POINTER(xout);
  youtP = NUMERIC_POINTER(yout);
  for(k = 0; k < nout; k++) {
    ioutP[k] = ii[k];
    joutP[k] = jj[k];
    tioutP[k] = ti[k];
    tjoutP[k] = tj[k];
    xoutP[k]  = x[k];
    youtP[k]  = y[k];
  }
  PROTECT(out = NEW_LIST(6));
  /* 11 + 1 = 12 protected */
  SET_VECTOR_ELT(out, 0, iout);
  SET_VECTOR_ELT(out, 1, jout);
  SET_VECTOR_ELT(out, 2, tiout);
  SET_VECTOR_ELT(out, 3, tjout);
  SET_VECTOR_ELT(out, 4, xout);
  SET_VECTOR_ELT(out, 5, yout);
  UNPROTECT(12);
  return(out);
}
