#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"


/* Conditional intensity computation for Soft Core process */

/* Storage of parameters and precomputed/auxiliary data */

typedef struct Softcore {
  double sigma;
  double kappa;
  double nook;  /*   -1/kappa     */
  double stok; /* sigma^(2/kappa) */
  double *period;
  int per;
} Softcore;


/* initialiser function */

Cdata *sftcrinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  Softcore *softcore;
  softcore = (Softcore *) R_alloc(1, sizeof(Softcore));

  /* Interpret model parameters*/
  softcore->sigma  = model.ipar[0];
  softcore->kappa  = model.ipar[1];
  softcore->period = model.period;
  /* constants */
  softcore->nook = -1/softcore->kappa;
  softcore->stok = pow(softcore->sigma, 2/softcore->kappa);
  /* periodic boundary conditions? */
  softcore->per    = (model.period[0] > 0.0);
  return((Cdata *) softcore);
}

/* conditional intensity evaluator */

double sftcrcif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, pairsum, cifval, nook, stok;
  Softcore *softcore;

  softcore = (Softcore *) cdata;

  nook = softcore->nook;
  stok = softcore->stok;

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
  if(softcore->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],softcore->period);
	pairsum += pow(d2, nook);
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],softcore->period);
	pairsum += pow(d2, nook);
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = pow(u - x[j],2) + pow(v-y[j],2);
	pairsum += pow(d2, nook);
      }
    }  
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = pow(u - x[j],2) + pow(v-y[j],2);
	pairsum += pow(d2, nook);
      }
    }
  }

  cifval *= exp(-stok * pairsum);
  return cifval;
}

Cifns SoftcoreCifns = { &sftcrinit, &sftcrcif, (updafunptr) NULL, FALSE};

