#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"
#include "constants.h"

/* Conditional intensity computation for Penttinen process */

/* Format for storage of parameters and precomputed/auxiliary data */

typedef struct Penttinen {
  double gamma;
  double r;
  double loggamma;
  double reach2;
  double *period;
  int hard;
  int per;
} Penttinen;


/* initialiser function */

Cdata *penttineninit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* create storage for model parameters */
  Penttinen *penttinen;
  penttinen = (Penttinen *) R_alloc(1, sizeof(Penttinen)); 
  /* Interpret model parameters*/
  penttinen->gamma  = model.ipar[0];
  penttinen->r      = model.ipar[1]; 
  penttinen->reach2 = 4.0 * penttinen->r * penttinen->r; 
  penttinen->period = model.period;
#ifdef MHDEBUG
  Rprintf("Initialising Penttinen gamma=%lf, r=%lf\n", 
	  penttinen->gamma, penttinen->r);
#endif
  /* is the model numerically equivalent to hard core ? */
  penttinen->hard   = (penttinen->gamma < DOUBLE_EPS);
  penttinen->loggamma = (penttinen->hard) ? 0 : log(penttinen->gamma);
  /* periodic boundary conditions? */
  penttinen->per    = (model.period[0] > 0.0);
  return((Cdata *) penttinen);
}

/* conditional intensity evaluator */

double penttinencif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, reach2, z, z2, logpot, cifval;
  Penttinen *penttinen;
  DECLARE_CLOSE_D2_VARS;

  penttinen = (Penttinen *) cdata;

  reach2     = penttinen->reach2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;

  if(npts == 0) 
    return((double) 1.0);

  logpot = 0.0;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(penttinen->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],penttinen->period,reach2,d2)) {
	  z2 = d2/reach2;
	  z = sqrt(z2);
	  if(z < 1.0) {
	    logpot += acos(z) - z * sqrt(1 - z2);
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],penttinen->period,reach2,d2)) {
	  z2 = d2/reach2;
	  z = sqrt(z2);
	  if(z < 1.0) {
	    logpot += acos(z) - z * sqrt(1 - z2);
	  }
	}
      }
    }
  } else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE_D2(u, v, x[j], y[j], reach2, d2)) {
	  z2 = d2/reach2;
	  z = sqrt(z2);
	  if(z < 1.0) {
	    logpot += acos(z) - z * sqrt(1 - z2);
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE_D2(u, v, x[j], y[j], reach2, d2)) {
	  z2 = d2/reach2;
	  z = sqrt(z2);
	  if(z < 1.0) {
	    logpot += acos(z) - z * sqrt(1 - z2);
	  }
	}
      }
    }
  }

  if(penttinen->hard) {
    if(logpot > 0) cifval = 0.0;
    else cifval = 1.0;
  } else cifval = exp((penttinen->loggamma) * M_2_PI * logpot);
  
  return cifval;
}

Cifns PenttinenCifns = { &penttineninit, &penttinencif, (updafunptr) NULL, NO};
