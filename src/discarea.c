/*

  disc.c

  Area of intersection between disc and polygonal window

  $Revision: 1.7 $     $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

#undef DEBUG

#include <math.h>

#include <R.h>

#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) > (B)) ? (A) : (B))

#ifndef PI
#define PI 3.1415926535898
#endif

void 
discareapoly(nc, xc, yc, nr, rmat, nseg, x0, y0, x1, y1, eps, out) 
     /* inputs */
     int *nc, *nr, *nseg;
     double *xc, *yc, *rmat;
     double *x0, *y0, *x1, *y1;
     double *eps;
     /* output */
     double *out;
{
  int n, m, i, j, k, nradperpt;
  double radius, radius2, total, contrib;
  double xx0, xx1, yy0, yy1, xleft, xright, yleft, yright, xcentre, ycentre;
  double epsilon;
  double DiscContrib();

  n = *nc;
  nradperpt = *nr;
  m = *nseg;
  epsilon = *eps;

  for(i = 0; i < n; i++) {
    xcentre = xc[i];
    ycentre = yc[i];
#ifdef DEBUG
    Rprintf("\ni = %d:\n centre = (%lf, %lf)\n", i, xcentre, ycentre);
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
	xx0 = x0[k];
	yy0 = y0[k];
	xx1 = x1[k];
	yy1 = y1[k];
#ifdef DEBUG
       Rprintf("(%lf,%lf) to (%lf,%lf)\n", xx0, yy0, xx1, yy1);
#endif
       /* refer to unit disc at origin */
       /* arrange so that xleft < xright */
       if(radius <= epsilon)
	 contrib = 0.0;
       else if(xx0 < xx1) {
	 xleft  = (xx0 - xcentre)/radius;
	 xright = (xx1 - xcentre)/radius;
	 yleft  = (yy0 - ycentre)/radius;
	 yright = (yy1 - ycentre)/radius;
	 contrib = - radius2 * DiscContrib(xleft,yleft,xright,yright,epsilon);
       } else {
	 xleft  = (xx1 - xcentre)/radius;
	 xright = (xx0 - xcentre)/radius;
	 yleft  = (yy1 - ycentre)/radius;
	 yright = (yy0 - ycentre)/radius;
	 contrib =  radius2 * DiscContrib(xleft,yleft,xright,yright,epsilon);
       }
#ifdef DEBUG
	Rprintf("contrib = %lf\n contrib/(pi * r^2)=%lf\n", 
		contrib, contrib/(PI * radius2));
#endif
	total += contrib;
      }
      out[ j * n + i] = total;
#ifdef DEBUG
	Rprintf("total = %lf\ntotal/(pi * r^2) = %lf\n", 
		total, total/(PI * radius2));
#endif
    }
  }
}

/* area of intersection of unit disc with halfplane x <= v */

#ifdef DEBUG
#define TRIGBIT(V) trigbit(V)
double trigbit(v) 
     double v;
{
  double zero, result;
  zero = 0.0;
  if(v < -1.0)
    return(zero);
  if(v > 1.0)
    return(PI);
  result = PI/2 + asin(v) + v * sqrt(1 - v * v);
  Rprintf("trigbit: v = %lf, asin(v)=%lf, result=%lf\n",
	  v, asin(v), result);
  return(result);
}
#else
#define TRIGBIT(V) (((V) <= -1.0) ? 0.0 : (((V) >= 1.0) ? PI : \
              (PI/2 + asin(V) + (V) * sqrt(1 - (V) * (V)))))
#endif

/*
  Find the area of intersection between a disc 
       centre = (0,0),   radius = 1
  and the trapezium with upper segment 
       (xleft, yleft) to (xright, yright)
  ASSUMES xleft < xright
*/

double DiscContrib(xleft, yleft, xright, yright, eps) 
  double xleft, yleft, xright, yright, eps;
  /* 
  NOTE: unit disc centred at origin
  */
{
  double xlo, xhi, zero, slope, intercept, A, B, C, det;
  double xcut1, xcut2, ycut1, ycut2, xunder1, xunder2, dx, dx2, result;

#ifdef DEBUG
  double increm;
  Rprintf(
	  "DiscContrib: xleft=%lf, yleft=%lf, xright=%lf, yright=%lf\n",
	  xleft, yleft, xright, yright);
#endif

  zero = 0.0;
  /* determine relevant range of x coordinates */
  xlo = MAX(xleft, (-1.0));
  xhi = MIN(xright, 1.0);
  if(xlo >= xhi - eps) {
    /* intersection is empty or negligible */
#ifdef DEBUG
    Rprintf("intersection is empty or negligible\n");
#endif
    return(zero);
  }
    
  /* find intersection points between the circle 
     and the line containing upper segment
  */
  slope = (yright - yleft)/(xright - xleft);
  intercept = yleft - slope * xleft;
  A = 1 + slope * slope;
  B = 2 * slope * intercept;
  C = intercept * intercept - 1.0;
  det = B * B - 4 * A * C;

#ifdef DEBUG
    Rprintf("slope=%lf, intercept=%lf\nA = %lf, B=%lf, C=%lf, det=%lf\n",
	    slope, intercept, A, B, C, det);
#endif

  if(det <= 0.0) {
    /* no crossing between disc and infinite line */
    if(intercept < 0.0) 
      /* segment is below disc; intersection is empty */
      return(zero);
    /* segment is above disc */
    result = TRIGBIT(xhi) - TRIGBIT(xlo);
    return(result);
  } 
  xcut1 = (- B - sqrt(det))/(2 * A);
  xcut2 = (- B + sqrt(det))/(2 * A);
  /* partition [xlo, xhi] into pieces delimited by {xcut1, xcut2} */
  if(xcut1 >= xhi || xcut2 <= xlo) {
    /* segment is outside disc */
    if(yleft < 0.0) {
#ifdef DEBUG
    Rprintf("segment is beneath disc\n");
#endif
      result = zero;
    } else {
#ifdef DEBUG
    Rprintf("segment is above disc\n");
#endif
      result = TRIGBIT(xhi) - TRIGBIT(xlo);
    }
    return(result);
  } 
  /* possibly three parts */
#ifdef DEBUG
  Rprintf("up to three pieces\n");
#endif
  result = zero;
  ycut1 = intercept + slope * xcut1;
  ycut2 = intercept + slope * xcut2;
  if(xcut1 > xlo) {
    /* part to left of cut */
#ifdef DEBUG 
    Rprintf("left of cut: [%lf, %lf]\n", xlo, xcut1);
    if(ycut1 < 0.0)
      Rprintf("below disc - no intersection\n");
    else {
      increm = TRIGBIT(xcut1) - TRIGBIT(xlo);
      Rprintf("increment = %lf\n", increm);
      result += increm;
    }
#else
    if(ycut1 >= 0.0)
      result += TRIGBIT(xcut1) - TRIGBIT(xlo);
#endif
  }
  if(xcut2 < xhi) {
    /* part to right of cut */
#ifdef DEBUG 
    Rprintf("right of cut: [%lf, %lf]\n", xcut2, xhi);
    if(ycut2 < 0.0)
      Rprintf("below disc - no intersection\n");
    else {
      increm = TRIGBIT(xhi) - TRIGBIT(xcut2);
      Rprintf("increment = %lf\n", increm);
      result += increm;
    }
#else
    if(ycut2 >= 0.0)
      result += TRIGBIT(xhi) - TRIGBIT(xcut2);
#endif
  }
  /* part underneath cut */
  xunder1 = MAX(xlo, xcut1);
  xunder2 = MIN(xhi, xcut2);
  dx = xunder2 - xunder1;
  dx2 = xunder2 * xunder2 - xunder1 * xunder1;
#ifdef DEBUG 
    Rprintf("underneath cut: [%lf, %lf]\n",
	    xunder1, xunder2);
    increm = intercept * dx + slope * dx2/2 + 
      (TRIGBIT(xunder2) - TRIGBIT(xunder1))/2;
    Rprintf("increment = %lf\n", increm);
    result += increm;
#else
  result += intercept * dx + slope * dx2/2 + 
    (TRIGBIT(xunder2) - TRIGBIT(xunder1))/2;
#endif
  
  return(result);

}
  

#ifdef DEBUG
/* interface to low level function, for debugging only */

void RDCtest(xleft, yleft, xright, yright, eps, value)
  double *xleft, *yleft, *xright, *yright, *eps, *value;
{
  double DiscContrib();
  *value = DiscContrib(*xleft, *yleft, *xright, *yright, *eps);
}

#endif
