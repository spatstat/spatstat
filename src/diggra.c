#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Diggle-Gratton process */

/*
 Conditional intensity function for a pairwise interaction point
 process with interaction function as given by 

                  e(t) = 0 for t < delta
                       = (t-delta)/(rho-delta)^kappa for delta <= t < rho
                       = 1 for t >= rho

 (See page 767 of Diggle, Gates, and Stibbard, Biometrika vol. 74,
  1987, pages 763 -- 770.)
*/

/* Storage of parameters and precomputed/auxiliary data */

typedef struct Diggra {
  double kappa;
  double delta;
  double rho;
  double delta2;  /*  delta^2   */
  double rho2;    /*  rho^2 */
  double fac;   /*   1/(rho-delta)  */
  double *period;
  int per;
} Diggra;


/* initialiser function */

Cdata *diggrainit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  Diggra *diggra;
  diggra = (Diggra *) R_alloc(1, sizeof(Diggra));

  /* Interpret model parameters*/
  diggra->kappa  = model.ipar[0];
  diggra->delta  = model.ipar[1];
  diggra->rho    = model.ipar[2];
  diggra->period = model.period;
  /* constants */
  diggra->delta2 = pow(diggra->delta, 2);
  diggra->rho2 = pow(diggra->rho, 2);
  diggra->fac = 1/(diggra->rho - diggra->delta);
  /* periodic boundary conditions? */
  diggra->per    = (model.period[0] > 0.0);
  return((Cdata *) diggra);
}

/* conditional intensity evaluator */

double diggracif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, pairprod, cifval;
  double rho2, delta, delta2, fac;
  double *period;
  DECLARE_CLOSE_D2_VARS;

  Diggra *diggra;

  diggra = (Diggra *) cdata;
  period = diggra->period;
  rho2   = diggra->rho2;
  delta  = diggra->delta;
  delta2 = diggra->delta2;
  fac    = diggra->fac;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = pairprod = 1.0;

  if(npts == 0) 
    return(cifval);

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(diggra->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,rho2,d2)) {
	  if(d2 < delta2) {
	    cifval = 0.0;
	    return(cifval);
	  } else {
	    pairprod *= fac * (sqrt(d2)-delta);
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,rho2,d2)) {
	  if(d2 < delta2) {
	    cifval = 0.0;
	    return(cifval);
	  } else {
	    pairprod *= fac * (sqrt(d2)-delta);
	  }
	}
      }
    }
  } else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
        if(CLOSE_D2(u, v, x[j], y[j], rho2, d2)) {
	  if(d2 <= delta2) {
	    cifval = 0.0;
	    return(cifval);
	  } else {
	    pairprod *= fac * (sqrt(d2)-delta);
	  }
	}
      }  
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
        if(CLOSE_D2(u, v, x[j], y[j], rho2, d2)) {
	  if(d2 <= delta2) {
	    cifval = 0.0;
	    return(cifval);
	  } else {
	    pairprod *= fac * (sqrt(d2)-delta);
	  }
	}
      }  
    }
  }

  cifval = pow(pairprod, diggra->kappa);
  return cifval;
}

Cifns DiggraCifns = { &diggrainit, &diggracif, (updafunptr) NULL, NO};

