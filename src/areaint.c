#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "methas.h"
#include "dist2.h"

/*
  Conditional intensity function for an area-interaction process:

  cif = eta^(1-B) where B = (uncovered area)/(pi r^2)

*/

#define NGRID 16

/* To explore serious bug, #define BADBUG */
#undef BADBUG

/* Format for storage of parameters and precomputed/auxiliary data */

typedef struct AreaInt {
  /* model parameters */
  double eta;
  double r;
  /* transformations of the parameters */
  double r2;
  double range2;
  double logeta;
  int hard;
  /* periodic distance */
  double *period;
  int per;
  /* grid counting */
  double dx;
  double xgrid0;
  int *my;
  int kdisc;
  /* scratch space for saving list of neighbours */
  int *neighbour;
} AreaInt;

/* initialiser function */

Cdata *areaintInit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  double r, dx, dy, x0;
  int i, my, kdisc;
  AreaInt *areaint;

  /* create storage */
  areaint = (AreaInt *) R_alloc(1, sizeof(AreaInt));
  /* Interpret model parameters*/
  areaint->eta    = model.ipar[0];
  areaint->r      = r = model.ipar[1]; 
#ifdef BADBUG
  Rprintf("r = %lf\n", r);
#endif
  areaint->r2     = r * r;
  areaint->range2 = 4 * r * r;    /* square of interaction distance */
  /* is the model numerically equivalent to hard core ? */
  areaint->hard   = (areaint->eta == 0.0);
  areaint->logeta = (areaint->hard) ? log(DOUBLE_XMIN) : log(areaint->eta);
#ifdef BADBUG
  if(areaint->hard) Rprintf("Hard core recognised\n");
#endif
  /* periodic boundary conditions? */
  areaint->period = model.period;
  areaint->per    = (model.period[0] > 0.0);
#ifdef BADBUG
  if(areaint->per) {
    Rprintf("*** periodic boundary conditions ***\n");
    Rprintf("period = %lf, %lf\n", model.period[0], model.period[1]);
  }
#endif
  /* grid counting */
  dx = dy = areaint->dx = (2 * r)/NGRID;
#ifdef BADBUG
  Rprintf("areaint->dx = %lf\n", areaint->dx);
#endif
  areaint->xgrid0 = -r + dx/2;
  areaint->my = (int *) R_alloc((long) NGRID, sizeof(int));
  kdisc = 0;
  for(i = 0; i < NGRID; i++) {
    x0 = areaint->xgrid0 + i * dx;
    my = floor(sqrt(r * r - x0 * x0)/dy);
    my = (my < 0) ? 0 : my;
    areaint->my[i] = my;
#ifdef BADBUG
    Rprintf("\tmy[%ld] = %ld\n", i, my);
#endif
    kdisc += 2 * my + 1;
  }
  areaint->kdisc = kdisc;
#ifdef BADBUG
  Rprintf("areaint->kdisc = %ld\n", areaint->kdisc);
#endif
  /* allocate space for neighbour indices */
  areaint->neighbour = (int *) R_alloc((long) state.npmax, sizeof(int));
  return((Cdata *) areaint);
}

#ifdef BADBUG
void fexitc();
#endif

/* conditional intensity evaluator */

double areaintCif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double r2, d2, dx, dy, a, range2;
  double xgrid, ygrid, xgrid0, ygrid0, covfrac, cifval;
  int kount, kdisc, kx, my, ky;
  int *neighbour;
  int nn, k;

  AreaInt *areaint;

  areaint = (AreaInt *) cdata;

  r2      = areaint->r2;
  range2  = areaint->range2;    /* square of interaction distance */
  dy = dx = areaint->dx;
  kdisc   = areaint->kdisc;
  /* pointers */
  period   = areaint->period;
  neighbour = areaint->neighbour;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  if(npts == 0) return ((double) 1.0);

  if(!areaint->per) {
    /*
      ..........   Euclidean distance ....................
      First identify which data points are neighbours of (u,v)
    */
    nn = 0;
    ixp1 = ix + 1;
    /* If ix = NONE = -1, then ixp1 = 0 is correct */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	a = range2 - pow(u - x[j], 2);
	if(a > 0.) {
	  a -= pow(v - y[j], 2);
	  if(a > 0.) {
	    /* point j is a neighbour of (u,v) */
	    neighbour[nn] = j;
	    ++nn;
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j < npts; j++) {
	a = range2 - pow(u - x[j], 2);
	if(a > 0.) {
	  a -= pow(v - y[j], 2);
	  if(a > 0.) {
	    /* point j is a neighbour of (u,v) */
	    neighbour[nn] = j;
	    ++nn;
	  }
	}
      }
    }
    if(nn == 0) {
      /* no neighbours; no interaction */
      cifval = 1.0;
      return cifval;
    } else if(areaint->hard) {
      /* neighbours forbidden if it's a hard core process */
      cifval = 0.0;
      return cifval;
    } else {
      /* scan a grid of points centred at (u,v) */
      kount = 0;
      xgrid0 = u + areaint->xgrid0;
      for(kx=0; kx<NGRID; kx++) {
	xgrid = xgrid0 + kx * dx;
	my = areaint->my[kx];
	for(ky=(-my); ky<=my; ky++) {
	  ygrid = v + ky * dy;
	  /*
	    Grid point (xgrid,ygrid) is inside disc of
	    radius r centred at (u,v)

	    Loop through all neighbouring data points to determine
	    whether the grid point is covered by another disc
	  */
	  if(nn > 0) {
	    for(k=0; k < nn; k++) {
	      j = neighbour[k];
	      a = r2 - pow(xgrid - x[j], 2);
	      if(a > 0) {
		a -= pow(ygrid - y[j], 2);
		if(a > 0) {
		  /* point j covers grid point */
		  ++kount;
		  break;
		}
	      }
	    }
	  }
	  /* finished consideration of grid point (xgrid, ygrid) */
	}
      }
    }
  } else {
    /*
      ............. periodic distance ......................
      First identify which data points are neighbours of (u,v)
    */
    nn = 0;
    ixp1 = ix + 1;
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(dist2thresh(u,v,x[j],y[j],period,range2)) {
	  /* point j is a neighbour of (u,v) */
	  neighbour[nn] = j;
	  ++nn;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(dist2thresh(u,v,x[j],y[j],period,range2)) {
	  /* point j is a neighbour of (u,v) */
	  neighbour[nn] = j;
	  ++nn;
	}
      }
    }
    if(nn == 0) {
      /* no neighbours; no interaction */
      cifval = 1.0;
      return cifval;
    } else if(areaint->hard) {
      /* neighbours forbidden if it's a hard core process */
      cifval = 0.0;
      return cifval;
    } else {
      /* scan a grid of points centred at (u,v) */
      kount = 0;
      xgrid0 = u + areaint->xgrid0;
      for(kx=0; kx<NGRID; kx++) {
	xgrid = xgrid0 + kx * dx;
	my = areaint->my[kx];
	for(ky=(-my); ky<=my; ky++) {
	  ygrid = v + ky * dy;
	  /*
	    Grid point (xgrid,ygrid) is inside disc of
	    radius r centred at (u,v)

	    Loop through all neighbouring data points to determine
	    whether the grid point is covered by another disc
	  */
	  for(k=0; k < nn; k++) {
	    j = neighbour[k];
	    if(dist2Mthresh(xgrid,ygrid,x[j],y[j],period,r2)) {  
	      /* point j covers grid point */
	      ++kount;
	      break;
	    }
	  }
	  /* finished considering grid point (xgrid,ygrid) */
	}
      }
    }
  }
  /*
    `kdisc' is the number of         grid points in the disc
    `kount' is the number of COVERED grid points in the disc
  */

  /* Hard core case has been handled. */
  /* Usual calculation: covered area fraction */
  covfrac = ((double) kount)/((double) kdisc);
  cifval = exp(areaint->logeta * covfrac);

#ifdef BADBUG
    if(!R_FINITE(cifval)) {
      Rprintf("Non-finite CIF value\n");
      Rprintf("kount=%ld, kdisc=%ld, covfrac=%lf, areaint->logeta=%lf\n", 
	      kount, kdisc, covfrac, areaint->logeta);
      Rprintf("u=%lf, v=%lf\n", u, v);
      fexitc("Non-finite CIF");
    }
#endif

  return cifval;
}


Cifns AreaIntCifns = { &areaintInit, &areaintCif, (updafunptr) NULL, FALSE};
