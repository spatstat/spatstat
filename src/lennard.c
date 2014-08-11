#include <R.h>
#include <Rmath.h>
#include <R_ext/Constants.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Lennard-Jones process */

/* Storage of parameters and precomputed/auxiliary data */

typedef struct Lennard {
  double sigma;
  double epsilon;
  double sigma2;  /*   sigma^2     */
  double foureps;    /*   4 * epsilon     */
  double d2min;  /* minimum value of d^2 which yields nonzero intensity */
  double d2max;  /* maximum value of d^2 which has nontrivial contribution */
  double *period;
  int per;
} Lennard;

/* 
   MAXEXP is intended to be the largest x such that exp(-x) != 0 
   although the exact value is not needed
*/
#define MAXEXP (-log(DOUBLE_XMIN))
#define MINEXP (log(1.001))

/* initialiser function */

Cdata *lennardinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  Lennard *lennard;
  double sigma2, foureps, minfrac, maxfrac;

  lennard = (Lennard *) R_alloc(1, sizeof(Lennard));

  /* Interpret model parameters*/
  lennard->sigma   = model.ipar[0];
  lennard->epsilon = model.ipar[1];
  lennard->period  = model.period;
  /* constants */
  lennard->sigma2  = sigma2 = pow(lennard->sigma, 2);
  lennard->foureps = foureps = 4 * lennard->epsilon;
  /* thresholds where the interaction becomes trivial */
  minfrac = pow(foureps/MAXEXP, (double) 1.0/6.0);
  if(minfrac > 0.5) minfrac = 0.5;
  maxfrac = pow(foureps/MINEXP, (double) 1.0/3.0);
  if(maxfrac < 2.0) maxfrac = 2.0;
  lennard->d2min   = sigma2 * minfrac;
  lennard->d2max   = sigma2 * maxfrac;
  /* periodic boundary conditions? */
  lennard->per    = (model.period[0] > 0.0);

  return((Cdata *) lennard);
}

/* conditional intensity evaluator */

double lennardcif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, ratio6, pairsum, cifval;
  double sigma2, d2max, d2min;
  double *period;
  Lennard *lennard;
  DECLARE_CLOSE_D2_VARS;

  lennard = (Lennard *) cdata;

  sigma2 = lennard->sigma2;
  d2max  = lennard->d2max;
  d2min  = lennard->d2min;
  period = lennard->period;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = 1.0;

  if(npts == 0) 
    return(cifval);

  pairsum = 0;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(lennard->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,d2max,d2)) {
	  if(d2 < d2min) {
	    cifval = 0.0;
	    return cifval;
	  }
	  ratio6 = pow(sigma2/d2, 3);
	  pairsum += ratio6 * (1.0 - ratio6);
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,d2max,d2)) {
	  if(d2 < d2min) {
	    cifval = 0.0;
	    return cifval;
	  }
	  ratio6 = pow(sigma2/d2, 3);
	  pairsum += ratio6 * (1.0 - ratio6);
	}
      }
    }
  } else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
        if(CLOSE_D2(u, v, x[j], y[j], d2max, d2)) {
	  if(d2 < lennard->d2min) {
	    cifval = 0.0;
	    return cifval;
	  }
	  ratio6 = pow(sigma2/d2, 3);
	  pairsum += ratio6 * (1.0 - ratio6);
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
        if(CLOSE_D2(u, v, x[j], y[j], d2max, d2)) {
	  if(d2 < lennard->d2min) {
	    cifval = 0.0;
	    return cifval;
	  }
	  ratio6 = pow(sigma2/d2, 3);
	  pairsum += ratio6 * (1.0 - ratio6);
	}
      }
    }
  }

  cifval *= exp(lennard->foureps * pairsum);
  return cifval;
}

Cifns LennardCifns = { &lennardinit, &lennardcif, (updafunptr) NULL, FALSE};

