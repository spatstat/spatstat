#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* for debugging code, include   #define DEBUG 1   */

/* Conditional intensity computation for Multitype Strauss hardcore process */

/* NOTE: types (marks) are numbered from 0 to ntypes-1 */

/* Storage of parameters and precomputed/auxiliary data */

typedef struct MultiStraussHard {
  int ntypes;
  double *gamma;   /* gamma[i,j] = gamma[i+ntypes*j] for i,j = 0... ntypes-1 */
  double *rad;     /* rad[i,j] = rad[j+ntypes*i] for i,j = 0... ntypes-1 */
  double *hc;      /* hc[i,j] = hc[j+ntypes*i] for i,j = 0... ntypes-1 */
  double *rad2;    /* squared radii */
  double *hc2;    /* squared radii */
  double *rad2hc2;    /* r^2 - h^2 */
  double  range2;   /* square of interaction range */
  double *loggamma; /* logs of gamma[i,j] */
  double *period;
  int    *hard;     /* hard[i,j] = 1 if gamma[i,j] ~~ 0 */
  int    *kount;    /* space for kounting pairs of each type */
  int per;
} MultiStraussHard;


/* initialiser function */

Cdata *straushminit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j, ntypes, n2, hard;
  double g, r, h, r2, h2, logg, range2;
  MultiStraussHard *multistrausshard;

  multistrausshard = (MultiStraussHard *) R_alloc(1, sizeof(MultiStraussHard));


  multistrausshard->ntypes = ntypes = model.ntypes;
  n2 = ntypes * ntypes;

#ifdef DEBUG
  Rprintf("initialising space for %d types\n", ntypes);
#endif

  /* Allocate space for parameters */
  multistrausshard->gamma    = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrausshard->rad      = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrausshard->hc       = (double *) R_alloc((size_t) n2, sizeof(double));

  /* Allocate space for transformed parameters */
  multistrausshard->rad2     = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrausshard->hc2      = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrausshard->rad2hc2  = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrausshard->loggamma = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrausshard->hard     = (int *) R_alloc((size_t) n2, sizeof(int));

  /* Allocate scratch space for counts of each pair of types */
  multistrausshard->kount    = (int *) R_alloc((size_t) n2, sizeof(int));

  /* Copy and process model parameters*/

  /* ipar will contain n^2 values of gamma, then n^2 values of r, 
     then n^2 values of h */

  range2 = 0.0;

  for(i = 0; i < ntypes; i++) {
    for(j = 0; j < ntypes; j++) {
      g = model.ipar[       i + j*ntypes];
      r = model.ipar[ n2  + i + j*ntypes];
      h = model.ipar[2*n2 + i + j*ntypes];
      r2 = r * r;
      h2 = h * h;
      hard = (g < DOUBLE_EPS);
      logg = (hard) ? 0 : log(g);
      MAT(multistrausshard->gamma,    i, j, ntypes) = g;
      MAT(multistrausshard->rad,      i, j, ntypes) = r;
      MAT(multistrausshard->hc,       i, j, ntypes) = h; 
      MAT(multistrausshard->rad2,     i, j, ntypes) = r2;
      MAT(multistrausshard->hc2,      i, j, ntypes) = h2;
      MAT(multistrausshard->rad2hc2,  i, j, ntypes) = r2-h2;
      MAT(multistrausshard->hard,     i, j, ntypes) = hard; 
      MAT(multistrausshard->loggamma, i, j, ntypes) = logg;
      if(r2 > range2) range2 = r2;
    }
  }
  multistrausshard->range2 = range2;

  /* periodic boundary conditions? */
  multistrausshard->period = model.period;
  multistrausshard->per    = (model.period[0] > 0.0);

#ifdef DEBUG
  Rprintf("end initialiser\n");
#endif
  return((Cdata *) multistrausshard);
}

/* conditional intensity evaluator */

double straushmcif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ntypes, kount, ix, ixp1, j, mrk, mrkj, m1, m2;
  int *marks;
  double *x, *y;
  double u, v, lg;
  double d2, a, cifval;
  double range2;
  double *period;
  MultiStraussHard *multistrausshard;
  DECLARE_CLOSE_D2_VARS;

  multistrausshard = (MultiStraussHard *) cdata;
  range2 = multistrausshard->range2;
  period = multistrausshard->period;

  u  = prop.u;
  v  = prop.v;
  mrk = prop.mrk;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  marks = state.marks;

  npts = state.npts;

#ifdef DEBUG
  Rprintf("computing cif: u=%lf, v=%lf, mrk=%d\n", u, v, mrk);
#endif

  cifval = 1.0;

  if(npts == 0) 
    return(cifval);

  ntypes = multistrausshard->ntypes;

#ifdef DEBUG
  Rprintf("initialising pair counts\n");
#endif

  /* initialise pair counts */
  for(m1 = 0; m1 < ntypes; m1++)
    for(m2 = 0; m2 < ntypes; m2++)
      MAT(multistrausshard->kount, m1, m2, ntypes) = 0;

  /* compile pair counts */

#ifdef DEBUG
  Rprintf("compiling pair counts\n");
#endif

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(multistrausshard->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,range2,d2)) {
	  mrkj = marks[j];
	  if(d2 < MAT(multistrausshard->rad2, mrk, mrkj, ntypes)) {
	    if(d2 < MAT(multistrausshard->hc2, mrk, mrkj, ntypes)) {
	      cifval = 0.0;
	      return(cifval);
	    }
	    MAT(multistrausshard->kount, mrk, mrkj, ntypes)++;
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	if(CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,range2,d2)) {
	  mrkj = marks[j];
	  if(d2 < MAT(multistrausshard->rad2, mrk, mrkj, ntypes)) {
	    if(d2 < MAT(multistrausshard->hc2, mrk, mrkj, ntypes)) {
	      cifval = 0.0;
	      return(cifval);
	    }
	    MAT(multistrausshard->kount, mrk, mrkj, ntypes)++;
	  }
	}
      }
    }
  } else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
        if(CLOSE_D2(u, v, x[j], y[j], range2, d2)) {
	  mrkj = marks[j];
	  if(d2 < MAT(multistrausshard->rad2, mrk, mrkj, ntypes)) {
	    if(d2 < MAT(multistrausshard->hc2, mrk, mrkj, ntypes)) {
	      cifval = 0.0;
	      return(cifval);
	    }
	    MAT(multistrausshard->kount, mrk, mrkj, ntypes)++;
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
        if(CLOSE_D2(u, v, x[j], y[j], range2, d2)) {
	  mrkj = marks[j];
	  if(d2 < MAT(multistrausshard->rad2, mrk, mrkj, ntypes)) {
	    if(d2 < MAT(multistrausshard->hc2, mrk, mrkj, ntypes)) {
	      cifval = 0.0;
	      return(cifval);
	    }
	    MAT(multistrausshard->kount, mrk, mrkj, ntypes)++;
	  }
	}
      }
    }
  }

#ifdef DEBUG
  Rprintf("multiplying cif factors\n");
#endif
  /* multiply cif value by pair potential */
  for(m1 = 0; m1 < ntypes; m1++) {
    for(m2 = 0; m2 < ntypes; m2++) {
      kount = MAT(multistrausshard->kount, m1, m2, ntypes);
      if(MAT(multistrausshard->hard, m1, m2, ntypes)) {
	if(kount > 0) {
	  cifval = 0.0;
	  return(cifval);
	}
      } else {
	lg = MAT(multistrausshard->loggamma, m1, m2, ntypes);
	cifval *= exp(lg * kount);
      }
    }
  }
  
#ifdef DEBUG
  Rprintf("returning positive cif\n");
#endif
  return cifval;
}

Cifns MultiStraussHardCifns = { &straushminit, &straushmcif, (updafunptr) NULL, YES};
