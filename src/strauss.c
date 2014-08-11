#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Strauss process */

/* Format for storage of parameters and precomputed/auxiliary data */

typedef struct Strauss {
  double gamma;
  double r;
  double loggamma;
  double r2;
  double *period;
  int hard;
  int per;
} Strauss;


/* initialiser function */

Cdata *straussinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* create storage for model parameters */
  Strauss *strauss;
  strauss = (Strauss *) R_alloc(1, sizeof(Strauss)); 
  /* Interpret model parameters*/
  strauss->gamma  = model.ipar[0];
  strauss->r      = model.ipar[1]; /* No longer passed as r^2 */
  strauss->r2     = strauss->r * strauss->r; 
  strauss->period = model.period;
#ifdef MHDEBUG
  Rprintf("Initialising Strauss gamma=%lf, r=%lf\n", 
	  strauss->gamma, strauss->r);
#endif
  /* is the model numerically equivalent to hard core ? */
  strauss->hard   = (strauss->gamma < DOUBLE_EPS);
  strauss->loggamma = (strauss->hard) ? 0 : log(strauss->gamma);
  /* periodic boundary conditions? */
  strauss->per    = (model.period[0] > 0.0);
  return((Cdata *) strauss);
}

/* conditional intensity evaluator */

double strausscif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, kount, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double r2, d2, a, cifval;
  Strauss *strauss;
  DECLARE_CLOSE_VARS;

  strauss = (Strauss *) cdata;

  r2     = strauss->r2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;

  if(npts == 0) 
    return((double) 1.0);

  kount = 0;
  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(strauss->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE_PERIODIC(u,v,x[j],y[j],strauss->period, r2))
	  ++kount;
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE_PERIODIC(u,v,x[j],y[j],strauss->period, r2))
	  ++kount;
      }
    }
  } else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE(u,v,x[j],y[j], r2))
	  ++kount;
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE(u,v,x[j],y[j], r2))
	  ++kount;
      }
    }
  }

  if(strauss->hard) {
    if(kount > 0) cifval = 0.0;
    else cifval = 1.0;
  }
  else cifval = exp((strauss->loggamma) * kount);
  
  return cifval;
}

Cifns StraussCifns = { &straussinit, &strausscif, (updafunptr) NULL, NO};
