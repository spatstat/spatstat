#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "methas.h"
#include "dist2.h"

void fexitc(const char *msg);

#undef MH_DEBUG 

/*
  Conditional intensity function for a Geyer saturation process.  
*/

typedef struct Geyer {
  /* model parameters */
  double gamma;
  double r;
  double s;
  /* transformations of the parameters */
  double r2;
  double loggamma;
  int hard;
  /* periodic distance */
  double *period;
  int per;
  /* auxiliary counts */
  int *aux;
#ifdef MH_DEBUG
  int *freshaux;
  int prevtype;
#endif
} Geyer;

Cdata *geyerinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j, n1;
  Geyer *geyer;
  double r2;
  double *period;
  DECLARE_CLOSE_VARS;

  geyer = (Geyer *) R_alloc(1, sizeof(Geyer));

  /* Interpret model parameters*/
  geyer->gamma  = model.ipar[0];
  geyer->r      = model.ipar[1]; /* not squared any more */
  geyer->s      = model.ipar[2]; 
  geyer->r2     = geyer->r * geyer->r;
#ifdef MHDEBUG
  Rprintf("Initialising Geyer gamma=%lf, r=%lf, sat=%lf\n",
	  geyer->gamma, geyer->r, geyer->s);
#endif
  /* is the model numerically equivalent to hard core ? */
  geyer->hard   = (geyer->gamma < DOUBLE_EPS);
  geyer->loggamma = (geyer->hard) ? 0 : log(geyer->gamma);
  /* periodic boundary conditions? */
  geyer->period = model.period;
  geyer->per    = (model.period[0] > 0.0);
  /* allocate storage for auxiliary counts */
  geyer->aux = (int *) R_alloc((size_t) state.npmax, sizeof(int));
#ifdef MH_DEBUG
  geyer->freshaux = (int *) R_alloc((size_t) state.npmax, sizeof(int));
  geyer->prevtype = -42;
#endif

  r2 = geyer->r2;

  /* Initialise auxiliary counts */
  for(i = 0; i < state.npmax; i++) 
    geyer->aux[i] = 0;

  if(geyer->per) {
    /* periodic */
    period = geyer->period;
    if(state.npts > 1) {
      n1 = state.npts - 1;
      for(i = 0; i < n1; i++) {
	for(j = i+1; j < state.npts; j++) {
	  if(CLOSE_PERIODIC(state.x[i], state.y[i], 
			    state.x[j], state.y[j], 
			    period, r2)) {
	    geyer->aux[i] += 1;
	    geyer->aux[j] += 1;
	  }
	}
      }
    }
  } else {
    /* Euclidean distance */
    if(state.npts > 1) {
      n1 = state.npts - 1;
      for(i = 0; i < n1; i++) {
	for(j = i+1; j < state.npts; j++) {
	  if(CLOSE(state.x[i], state.y[i], 
		 state.x[j], state.y[j], 
		   r2)) {
	    geyer->aux[i] += 1;
	    geyer->aux[j] += 1;
	  }
	}
      }
    }
  }
  return((Cdata *) geyer);
}

double geyercif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int ix, j, npts, tee;
  double u, v, r2, s;
  double w, a, b, f, cifval;
  double *x, *y;
  int *aux;
  double *period;
  Geyer *geyer;
  DECLARE_CLOSE_VARS;

  geyer = (Geyer *) cdata;

  npts = state.npts;
  if(npts==0) return ((double) 1.0);

  x = state.x;
  y = state.y;
  u = prop.u;
  v = prop.v;
  ix = prop.ix;

  r2     = geyer->r2;
  s      = geyer->s;
  period = geyer->period;
  aux    = geyer->aux;

  /* 
     tee = neighbour count at the point in question;
     w   = sum of changes in (saturated) neighbour counts at other points 
  */
  tee = w = 0.0;

  if(prop.itype == BIRTH) {
    if(geyer->per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	if(CLOSE_PERIODIC(u,v,x[j],y[j],period,r2)) {
	  tee++;
	  f = s - aux[j];
	  if(f > 1) /* j is not saturated after addition of (u,v) */
	    w = w + 1; /* addition of (u,v) increases count by 1 */
	  else if(f > 0) /* j becomes saturated by addition of (u,v) */
	    w = w + f;
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(CLOSE(u,v,x[j],y[j],r2)) {
	  tee++;
	  f = s - aux[j];
	  if(f > 1) /* j is not saturated after addition of (u,v) */
	    w = w + 1; /* addition of (u,v) increases count by 1 */
	  else if(f > 0) /* j becomes saturated by addition of (u,v) */
	    w = w + f;
	}
      }
    }
  } else if(prop.itype == DEATH) {
    tee = aux[ix];
    if(geyer->per) {
      /* Periodic distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	if(CLOSE_PERIODIC(u,v,x[j],y[j],period,r2)) {
	  f = s - aux[j];
	  if(f > 0) /* j is not saturated */
	    w = w + 1; /* deletion of 'ix' decreases count by 1 */
	  else {
	    f = f+1;
	    if(f > 0) {
	      /* j is not saturated after deletion of 'ix' 
		 (s must be fractional) */
	      w = w + f; 
	    }
	  }
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	if(CLOSE(u,v,x[j],y[j],r2)) {
	  f = s - aux[j];
	  if(f > 0) /* j was not saturated */
	    w = w + 1; /* deletion of 'ix' decreases count by 1 */
	  else {
	    f = f+1; 
	    if(f > 0) {
	      /* j is not saturated after deletion of 'ix' 
		 (s must be fractional) */
	      w = w + f; 
	    }
	  }
	}
      }
    }
  } else if(prop.itype == SHIFT) { 
    /* Compute the cif at the new point, not the ratio of new/old */
    if(geyer->per) {
      /* Periodic distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	if(CLOSE_PERIODIC(u,v,x[j],y[j],period,r2)) {
	  tee++;
	  a = aux[j];
	  /* Adjust */
	  if(CLOSE_PERIODIC(x[ix],y[ix],x[j],y[j],period,r2)) a = a - 1;
	  b = a + 1;
	  if(a < s && s < b) {
	    w = w + s - a;
	  }
	  else if(s >= b) w = w + 1;
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	if(CLOSE(u,v,x[j],y[j],r2)) {
	  tee++;
	  a = aux[j];
	  /* Adjust */
	  if(CLOSE(x[ix], y[ix], x[j], y[j], r2)) a = a - 1;
	  b = a + 1;
	  if(a < s && s < b) {
	    w = w + s - a;
	  }
	  else if(s >= b) w = w + 1;
	}
      }
    }
  }

  w = w + ((tee < s) ? tee : s);

 if(geyer->hard) {
    if(tee > 0) cifval = 0.0;
    else cifval = 1.0;
  }
  else cifval = exp(geyer->loggamma*w);
  
  return cifval;
}

void geyerupd(state, prop, cdata) 
     State state;
     Propo prop;
     Cdata *cdata;
{
/* Declare other variables */
  int ix, npts, j;
  int oldclose, newclose;
  double u, v, xix, yix, r2;
  double *x, *y;
  int *aux;
  double *period;
  Geyer *geyer;
#ifdef MH_DEBUG
  int *freshaux;
  int i;
  int oc, nc;
#endif
  DECLARE_CLOSE_VARS;

  geyer = (Geyer *) cdata;
  period = geyer->period;
  aux = geyer->aux;
  r2 = geyer->r2;

  x = state.x;
  y = state.y;
  npts = state.npts;

#ifdef MH_DEBUG  
  /* ........................ debugging cross-check ................ */

  /* recompute 'aux' values afresh */
  freshaux = geyer->freshaux;
  for(i = 0; i < state.npts; i++)
    freshaux[i] = 0;

  if(geyer->per) {
    /* periodic */
    for(i = 0; i < state.npts; i++) {
      for(j = 0; j < state.npts; j++) {
	if(i == j) continue;
	if(CLOSE_PERIODIC(state.x[i], state.y[i],
			  state.x[j], state.y[j],
			  period, r2)) 
	  freshaux[i] += 1;
      }
    }
  } else {
    /* Euclidean distance */
    for(i = 0; i < state.npts; i++) {
      for(j = 0; j < state.npts; j++) {
	if(i == j) continue;
	if(CLOSE(state.x[i], state.y[i], 
		 state.x[j], state.y[j], 
		 r2))
	  freshaux[i] += 1;
      }
    }
  }
  /* Check agreement with 'aux' */
  for(j = 0; j < state.npts; j++) {
    if(aux[j] != freshaux[j]) {
      Rprintf("\n\taux[%d] = %d, freshaux[%d] = %d\n", 
	      j, aux[j], j, freshaux[j]);
      Rprintf("\tnpts = %d\n", state.npts);
      Rprintf("\tperiod = (%lf, %lf)\n", period[0], period[1]);
      if(geyer->prevtype == BIRTH) error("updaux failed after BIRTH");
      if(geyer->prevtype == DEATH) error("updaux failed after DEATH");
      if(geyer->prevtype == SHIFT) error("updaux failed after SHIFT");
      error("updaux failed at start");
    }
  }
  /* OK. Record type of this transition */ 
  geyer->prevtype = prop.itype;

  /* ................ end debug cross-check ................ */
#endif

  if(prop.itype == BIRTH) { 
    /* Birth */
    u = prop.u;
    v = prop.v;
    /* initialise auxiliary counter for new point */
    aux[npts] = 0; 
    /* update all auxiliary counters */
    if(geyer->per) {
      /* periodic distance */
      for(j=0; j < npts; j++) {
	if(CLOSE_PERIODIC(u,v,x[j],y[j],period,r2)) {
	  aux[j] += 1;
	  aux[npts] += 1;
	} 
      }
    } else {
      /* Euclidean distance */
      for(j=0; j < npts; j++) {
	if(CLOSE(u,v,x[j],y[j],r2)) {
	  aux[j] += 1;
	  aux[npts] += 1;
	} 
      }
    }
  } else if(prop.itype == DEATH) {
    /* Death */
    ix = prop.ix;
    u = x[ix];
    v = y[ix];
    /* decrement auxiliary counter for each point */
    if(geyer->per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	if(j==ix) continue;
	if(CLOSE_PERIODIC(u,v,x[j],y[j],period,r2)) {
	  if(j < ix) aux[j] -= 1;
	  else aux[j-1] = aux[j] - 1;
	} else if(j >= ix) aux[j-1] = aux[j];
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j==ix) continue;
	if(CLOSE(u,v,x[j],y[j],r2)) {
	  if(j < ix) aux[j] -= 1;
	  else aux[j-1] = aux[j] - 1;
	} else if(j >= ix) aux[j-1] = aux[j];
      }
    }
  } else if(prop.itype == SHIFT) { 
    /* Shift */
    u = prop.u;
    v = prop.v;
    ix = prop.ix;
    xix = x[ix];
    yix = y[ix];
    /* recompute auxiliary counter for point 'ix' */
    aux[ix] = 0;
    /* update auxiliary counters for other points */
    if(geyer->per) {
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	newclose = oldclose = NO;
	if(CLOSE_PERIODIC(u,v,x[j],y[j],period,r2)) newclose = YES;
	if(CLOSE_PERIODIC(xix,yix,x[j],y[j],period,r2)) oldclose = YES;
	if(newclose) {
	  /* increment neighbour count for new point */
	  aux[ix] += 1;
	  if(!oldclose) 
	    aux[j] += 1; /* point j gains a new neighbour */
	} else if(oldclose)
	  aux[j] -= 1; /* point j loses a neighbour */
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	newclose = oldclose = NO;
	if(CLOSE(u,v,x[j],y[j],r2)) newclose = YES;
	if(CLOSE(xix,yix,x[j],y[j],r2)) oldclose = YES;
	if(newclose) {
	  /* increment neighbour count for new point */
	  aux[ix] += 1;
	  if(!oldclose) 
	    aux[j] += 1; /* point j gains a new neighbour */
	} else if(oldclose)
	  aux[j] -= 1; /* point j loses a neighbour */
      }
    }
  } else fexitc("Unrecognised transition type; bailing out.\n");

  return;
}

Cifns GeyerCifns = { &geyerinit, &geyercif, &geyerupd, NO};
