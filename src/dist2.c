# include <math.h>
#include <R.h>

#include "yesno.h"

/* 

   dist2:   squared distance in torus

   dist2thresh: faster code for testing whether dist2 < r2

   dist2Mthresh: same as dist2thresh, but does not assume
                 the points are within one period of each other.

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

double dist2(u,v,x,y,period)
     double u, v, x, y;
     double *period;
{
  double wide, high, dx, dy, dxp, dyp, a, b, d2;
  /* points are assumed to lie within one period of each other */

  wide = period[0];
  high = period[1];

  dx = u - x;
  if(dx < 0.0) dx = -dx;
  dxp = wide - dx;
  a = (dx < dxp)? dx : dxp;

  dy = v - y;
  if(dy < 0.0) dy = -dy;
  dyp = high - dy;
  b = (dy < dyp)? dy : dyp;

  d2 = a * a + b * b;
  return d2;
}

double dist2either(u,v,x,y,period)
     double u, v, x, y;
     double *period;
{
  if(period[0] < 0.0) return pow(u-x,2) + pow(v-y,2);
  return(dist2(u,v,x,y,period));
}

int dist2thresh(u,v,x,y,period,r2)
     double u, v, x, y, r2;
     double *period;
{
  double wide, high, dx, dy, dxp, dyp, a, b, residue;
  /* points are assumed to lie within one period of each other */

  wide = period[0];
  high = period[1];

  dx = u - x;
  if(dx < 0.0) dx = -dx;
  dxp = wide - dx;
  a = (dx < dxp) ? dx : dxp;
  residue = r2 - a * a;
  if(residue <= 0.0)
    return NO;
  dy = v - y;
  if(dy < 0.0) dy = -dy;
  dyp = high - dy;
  b = (dy < dyp) ? dy : dyp;
  if (residue > b * b) 
    return YES; 
  return NO;
}

int dist2Mthresh(u,v,x,y,period,r2)
     double u, v, x, y, r2;
     double *period;
{
  double wide, high, dx, dy, dxp, dyp, a, b, residue;
  /* points are NOT assumed to lie within one period of each other */

  wide = period[0];
  high = period[1];

  dx = u - x;
  if(dx < 0.0) dx = -dx;
  while(dx > wide) dx -= wide;
  dxp = wide - dx;
  a = (dx < dxp) ? dx : dxp;
  residue = r2 - a * a;
  if(residue < 0.0)
    return NO;
  dy = v - y;
  if(dy < 0.0) dy = -dy;
  while(dy > high) dy -= high;
  dyp = high - dy;
  b = (dy < dyp) ? dy : dyp;
  if (residue >= b * b) 
    return YES; 
  return NO;
}
