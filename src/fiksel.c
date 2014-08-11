#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Fiksel process */

/*
 Conditional intensity function for a pairwise interaction point
 process with interaction function 

                  e(t) = 0 for t < h
                       = exp(a * exp(- kappa * t)) for h <= t < r
                       = 1 for t >= r

*/

/* Storage of parameters and precomputed/auxiliary data */

typedef struct Fiksel {
  double r;
  double h;
  double kappa;
  double a;
  double h2;  /*  h^2   */
  double r2;  /*  r^2 */
  double *period;
  int per;
} Fiksel;


/* initialiser function */

Cdata *fikselinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  Fiksel *fiksel;
  fiksel = (Fiksel *) R_alloc(1, sizeof(Fiksel));

  /* Interpret model parameters*/
  fiksel->r      = model.ipar[0];
  fiksel->h      = model.ipar[1];
  fiksel->kappa  = model.ipar[2];
  fiksel->a      = model.ipar[3];
  fiksel->period = model.period;
  /* constants */
  fiksel->h2 = pow(fiksel->h, 2);
  fiksel->r2 = pow(fiksel->r, 2);
  /* periodic boundary conditions? */
  fiksel->per    = (model.period[0] > 0.0);

  return((Cdata *) fiksel);
}

/* conditional intensity evaluator */

double fikselcif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, pairpotsum, cifval;
  double kappa, r2, h2;
  double *period;
  Fiksel *fiksel;
  DECLARE_CLOSE_D2_VARS;

  fiksel = (Fiksel *) cdata;
  period = fiksel->period;
  kappa  = fiksel->kappa;
  r2     = fiksel->r2;
  h2     = fiksel->h2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = 1.0;

  if(npts == 0) 
    return(cifval);

  pairpotsum = 0;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(fiksel->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,r2,d2)) {	
	  if(d2 < h2) {
	    cifval = 0.0;
	    return(cifval);
	  } else {
	    pairpotsum += exp(-kappa * sqrt(d2));
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,r2,d2)) {	
	  if(d2 < h2) {
	    cifval = 0.0;
	    return(cifval);
	  } else {
	    pairpotsum += exp(-kappa * sqrt(d2));
	  }
	}
      }
    }
  } else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE_D2(u,v,x[j],y[j],r2,d2)) {	
	  if(d2 < h2) {
	    cifval = 0.0;
	    return(cifval);
	  } else {
	    pairpotsum += exp(-kappa * sqrt(d2));
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE_D2(u,v,x[j],y[j],r2,d2)) {	
	  if(d2 < h2) {
	    cifval = 0.0;
	    return(cifval);
	  } else {
	    pairpotsum += exp(-kappa * sqrt(d2));
	  }
	}
      }
    }
  }

  cifval = exp(fiksel->a * pairpotsum);
  return cifval;
}

Cifns FikselCifns = { &fikselinit, &fikselcif, (updafunptr) NULL, FALSE};

