#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "methas.h"
#include "dist2.h"

/* To get debug output, insert the line:  #define DEBUG 1  */

void fexitc(const char *msg);

/*
  Conditional intensity function for a multiscale saturation process. 

  parameter vector: 
      ipar[0] = ndisc
      ipar[1] = gamma[0]
      ipar[2] = r[0]
      ipar[3] = s[0]
      ...
*/

typedef struct BadGey {
  /* model parameters */
  int ndisc;
  double *gamma;
  double *r;
  double *s;
  /* transformations of the parameters */
  double *r2;
  double *loggamma;
  int *hard;
  /* periodic distance */
  double *period;
  int per;
  /* auxiliary counts */
  int *aux;   /* matrix[ndisc, npmax]: neighbour counts in current state */
  int *tee;   /* vector[ndisc] : neighbour count at point in question */
  double *w;  /* vector[ndisc] : sum of changes in counts at other points */
} BadGey;

Cdata *badgeyinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j, k, i0, ndisc, nmatrix;
  double r, g, d2;
  BadGey *badgey;

  /* create storage */
  badgey = (BadGey *) R_alloc(1, sizeof(BadGey));

  badgey->ndisc  = ndisc = model.ipar[0];
  /* Allocate space for parameter vectors */
  badgey->gamma    = (double *) R_alloc((size_t) ndisc, sizeof(double));
  badgey->r        = (double *) R_alloc((size_t) ndisc, sizeof(double));
  badgey->s        = (double *) R_alloc((size_t) ndisc, sizeof(double));
  /* Derived values */
  badgey->r2       = (double *) R_alloc((size_t) ndisc, sizeof(double));
  badgey->loggamma = (double *) R_alloc((size_t) ndisc, sizeof(double));
  badgey->hard     = (int *) R_alloc((size_t) ndisc, sizeof(int));
  /* copy and transform parameters */
  for(i=0; i < ndisc; i++) {
    i0 = 3*i + 1;
    g = badgey->gamma[i] = model.ipar[i0];
    r = badgey->r[i] =     model.ipar[i0 + 1];
        badgey->s[i] =     model.ipar[i0 + 2];
    badgey->r2[i] = r * r;
    badgey->hard[i] = (g < DOUBLE_EPS);
    badgey->loggamma[i] = (g < DOUBLE_EPS) ? 0 : log(g);
  }
  /* periodic boundary conditions? */
  badgey->period = model.period;
  badgey->per    = (model.period[0] > 0.0);
  /* Allocate scratch space */
  badgey->tee      = (int *) R_alloc((size_t) ndisc, sizeof(int));
  badgey->w        = (double *) R_alloc((size_t) ndisc, sizeof(double));
  /* Allocate space for auxiliary counts */
  nmatrix = ndisc * state.npmax;
  badgey->aux      = (int *) R_alloc((size_t) nmatrix, sizeof(int));
  /* Initialise auxiliary counts */
  for(i = 0; i < nmatrix; i++)
    badgey->aux[i] = 0;
  for(i = 0; i < state.npts; i++) {
    for(j = 0; j < state.npts; j++) {
      if(j == i) continue;
      d2 = dist2either(state.x[i], state.y[i], state.x[j], state.y[j], 
		       badgey->period);
      for(k = 0; k < ndisc; k++) {
	if(d2 < badgey->r2[k])
	  MAT(badgey->aux, k, i, ndisc) += 1;
      }
    }
  }
#ifdef DEBUG
  Rprintf("Finished initialiser; ndisc=%d\n", ndisc);
#endif
  return((Cdata *) badgey);
}

#define AUX(I,J) MAT(aux, I, J, ndisc)

double badgeycif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int ix, j, k, npts, ndisc, tk;
  double u, v, d2;
  double a, dd2, b, f, r2, s, cifval;
  double *x, *y;
  int *tee, *aux;
  double *w;
  BadGey *badgey;

  badgey = (BadGey *) cdata;

#ifdef DEBUG
  Rprintf("Entering badgeycif\n");
#endif

  npts = state.npts;
  cifval = 1.0;
  if(npts==0) return cifval;

  x = state.x;
  y = state.y;
  u = prop.u;
  v = prop.v;
  ix = prop.ix;

  ndisc = badgey->ndisc;
  tee   = badgey->tee;
  aux   = badgey->aux;
  w     = badgey->w;

  /* 
     For disc k, 
     tee[k] = neighbour count at the point in question;
     w[k]   = sum of changes in (saturated) neighbour counts at other points 
  */
  if(prop.itype == BIRTH) {
    /* compute tee[k] and w[k] from scratch */
    for(k = 0; k < ndisc; k++) {
      tee[k] = 0;
      w[k] = 0.0;
    }
    if(badgey->per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],badgey->period);
	for(k = 0; k < ndisc; k++) {
	  if(d2 < badgey->r2[k]) {
	    tee[k]++;
	    f = badgey->s[k] - AUX(k,j);
	    if(f > 1) /* j is not saturated after addition of (u,v) */
	      w[k] += 1; /* addition of (u,v) increases count by 1 */
	    else if(f > 0) /* j becomes saturated by addition of (u,v) */
	      w[k] += f;
	  }
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	for(k = 0; k < ndisc; k++) {
	  if(d2 < badgey->r2[k]) {
	    tee[k]++;
	    f = badgey->s[k] - AUX(k,j);
	    if(f > 1) /* j is not saturated after addition of (u,v) */
	      w[k] += 1; /* addition of (u,v) increases count by 1 */
	    else if(f > 0) /* j becomes saturated by addition of (u,v) */
	      w[k] += f;
	  }
	}
      }
    }
  } else if(prop.itype == DEATH) {
    /* extract current auxiliary counts for point ix */
    /* compute w[k] from scratch */
    for(k = 0; k < ndisc; k++) {
      tee[k] = AUX(k,ix);
      w[k] = 0.0;
    }
    /* compute change in counts for other points */
    if(badgey->per) {
      /* Periodic distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2 = dist2(u,v,x[j],y[j],badgey->period);
	for(k = 0; k < ndisc; k++) {
	  if(d2 < badgey->r2[k]) {
	    f = badgey->s[k] - AUX(k,j);
	    if(f > 0) /* j is not saturated */
	      w[k] += 1; /* deletion of 'ix' decreases count by 1 */
	    else {
	      f += 1;
	      if(f > 0) {
		/* j is not saturated after deletion of 'ix' 
		   (s must be fractional) */
		w[k] += f; 
	      }
	    }
	  }
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	for(k = 0; k < ndisc; k++) {
	  if(d2 < badgey->r2[k]) {
	    f = badgey->s[k] - AUX(k,j);
	    if(f > 0) /* j is not saturated */
	      w[k] += 1; /* deletion of 'ix' decreases count by 1 */
	    else {
	      f += 1;
	      if(f > 0) {
		/* j is not saturated after deletion of 'ix' 
		   (s must be fractional) */
		w[k] += f; 
	      }
	    }
	  }
	}
      }
    }
  } else if(prop.itype == SHIFT) { 
    /* compute auxiliary counts from scratch */
    for(k = 0; k < ndisc; k++) {
      tee[k] = 0;
      w[k] = 0.0;
    }
    /* Compute the cif at the new point, not the ratio of new/old */
    if(badgey->per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2 = dist2(u,v,x[j],y[j],badgey->period);
	for(k = 0; k < ndisc; k++) {
	  r2 = badgey->r2[k];
	  if(d2 < r2) {
	    /* shifted point is a neighbour of point j */
	    tee[k]++;
	    a = AUX(k,j);
	    s = badgey->s[k];
	    /* Adjust */
	    dd2 = dist2(x[ix],y[ix], x[j],y[j],badgey->period);
	    if(dd2 < r2) a -= 1; 
	    b = a + 1;
	    /* b is the number of neighbours of point j in new state */
	    if(a < s && s < b) {
	      w[k] += s - a;  /* s is fractional and j is saturated */
	    }
	    else if(s >= b) w[k] += 1; 
	  }
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	for(k = 0; k < ndisc; k++) {
	  r2 = badgey->r2[k];
	  if(d2 < r2) {
	    /* shifted point is a neighbour of point j */
	    tee[k]++;
	    a = AUX(k,j);
	    s = badgey->s[k];
	    /* Adjust */
	    dd2 = pow(x[ix] - x[j], 2) + pow(y[ix] - y[j], 2);
	    if(dd2 < r2) a -= 1; 
	    b = a + 1;
	    /* b is the number of neighbours of point j in new state */
	    if(a < s && s < b) {
	      w[k] += s - a;  /* s is fractional and j is saturated */
	    }
	    else if(s >= b) w[k] += 1; 
	  }
	}
      }
    }
  }

#ifdef DEBUG
  Rprintf("ndisc=%d\n", ndisc);
#endif

  /* compute total change in saturated count */
  for(k = 0; k < ndisc; k++) {
    s = badgey->s[k];
    tk = tee[k];
    w[k] += ((tk < s) ? tk : s);
#ifdef DEBUG
    Rprintf("s[%d]=%lf, t[%d]=%d, w[%d]=%lf\n",
	   k, s, k, tk, k, w[k]);
#endif
  }

  /* evaluate cif */
  for(k = 0; k < ndisc; k++) {
    if(badgey->hard[k]) {
      if(tee[k] > 0) return(0.0);
      /* else cifval multiplied by 0^0 = 1 */
    } else cifval *= exp(badgey->loggamma[k] * w[k]);
  }
  
  return cifval;
}

void badgeyupd(state, prop, cdata) 
     State state;
     Propo prop;
     Cdata *cdata;
{
/* Declare other variables */
  int ix, npts, ndisc, j, k;
  double u, v, xix, yix, r2, d2, d2old, d2new;
  double *x, *y;
  int *aux;
  BadGey *badgey;

  badgey = (BadGey *) cdata;

  aux = badgey->aux;
  /* 'state' is current state before transition */
  x = state.x;
  y = state.y;
  npts = state.npts;      
  ndisc = badgey->ndisc;

#ifdef DEBUG
  Rprintf("start update ---- \n");
  for(j=0; j < npts; j++) {
    for(k=0; k < ndisc; k++)
      Rprintf("aux[%d,%d]=%d\t", k, j, AUX(k,j));
    Rprintf("\n");
  }
#endif
      
  if(prop.itype == BIRTH) { 
#ifdef DEBUG
    Rprintf("Update for birth ---- \n");
#endif
    /* Birth */
    u = prop.u;
    v = prop.v;
    /* initialise auxiliary counters for new point x[npts], y[npts] */
    for(k = 0; k < ndisc; k++)
      AUX(k, npts) = 0;
    /* update all auxiliary counters */
    if(badgey->per) {
      /* periodic distance */
      for(j=0; j < npts; j++) {
	d2 = dist2(u,v,x[j],y[j],badgey->period);
	for(k = 0; k < ndisc; k++) {
	  if(d2 < badgey->r2[k]) {
	    AUX(k, j) += 1;
	    AUX(k, npts) += 1;
	  }
	} 
      }
    } else {
      /* Euclidean distance */
      for(j=0; j < npts; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	for(k = 0; k < ndisc; k++) {
	  if(d2 < badgey->r2[k]) {
	    AUX( k, j) += 1;
	    AUX( k, npts) += 1;
	  }
	} 
      }
    }
#ifdef DEBUG
  Rprintf("end update ---- \n");
  for(j=0; j <= npts; j++) {
    for(k=0; k < ndisc; k++)
      Rprintf("aux[%d,%d]=%d\t", k, j, AUX(k,j));
    Rprintf("\n");
  }
#endif
    return;
  }
  if(prop.itype == DEATH) {
    /* Death */
    ix = prop.ix;
    u = x[ix];
    v = y[ix];
#ifdef DEBUG
    Rprintf("--- Update for death of point %d = (%lf,%lf) ---- \n", ix, u, v);
#endif
    /* 
       Decrement auxiliary counter for each neighbour of deleted point,
       and remove entry corresponding to deleted point
    */
    if(badgey->per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	if(j==ix) continue;
	d2 = dist2(u,v,x[j],y[j],badgey->period);
	for(k = 0; k < ndisc; k++) {
	  if(d2 < badgey->r2[k]) {
	    if(j < ix) AUX(k,j) -= 1; 
	    else AUX(k,j-1) = AUX(k,j) - 1;
	  } else if(j >= ix) AUX(k,j-1) = AUX(k,j);
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j==ix) continue;
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	for(k = 0; k < ndisc; k++) {
	  if(d2 < badgey->r2[k]) {
#ifdef DEBUG
	    Rprintf("hit for point %d with radius r[%d]\n", j, k);
#endif
	    if(j < ix) AUX(k,j) -= 1; 
	    else AUX(k,j-1) = AUX(k,j) - 1;
	  } else if(j >= ix) AUX(k,j-1) = AUX(k,j);
	}
      }
    }
#ifdef DEBUG
  Rprintf("end update ---- \n");
  for(j=0; j < npts-1; j++) {
    for(k=0; k < ndisc; k++)
      Rprintf("aux[%d,%d]=%d\t", k, j, AUX(k,j));
    Rprintf("\n");
  }
#endif
    return;
  }

  if(prop.itype == SHIFT) { 
#ifdef DEBUG
    Rprintf("Update for shift ---- \n");
#endif
    /* Shift */
    u = prop.u;
    v = prop.v;
    ix = prop.ix;
    xix = x[ix];
    yix = y[ix];
    /* recompute all auxiliary counters for point ix */
    for(k = 0; k < ndisc; k++) 
      AUX(k,ix) = 0;

    if(badgey->per) {
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2new = dist2(u,v,x[j],y[j],badgey->period);
	d2old = dist2(xix,yix,x[j],y[j],badgey->period);
	for(k = 0; k < ndisc; k++) {
	  r2 = badgey->r2[k];
	  if(d2old >= r2 && d2new >= r2) continue;
	  if(d2new < r2) {
	    /* increment neighbour count for new point */
	    AUX(k,ix) += 1;
	    if(d2old >= r2) 
	      AUX(k,j) += 1; /* point j gains a new neighbour */
	  } else if(d2old < r2) 
	    AUX(k,j) -= 1; /* point j loses a neighbour */
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2new = pow(u - x[j], 2) + pow(v - y[j], 2);
	d2old = pow(x[ix] - x[j], 2) + pow(y[ix] - y[j], 2);
	for(k = 0; k < ndisc; k++) {
	  r2 = badgey->r2[k];
	  if(d2old >= r2 && d2new >= r2) continue;
	  if(d2new < r2) {
#ifdef DEBUG
	    Rprintf("shifted point is close to j=%d\n", j);
#endif
	    /* increment neighbour count for new point */
	    AUX(k,ix) += 1;
	    if(d2old >= r2) {
#ifdef DEBUG
	    Rprintf("\t(previous position was not)\n");
#endif
	      AUX(k,j) += 1; /* point j gains a new neighbour */
	    }
	  } else if(d2old < r2) {
#ifdef DEBUG
	    Rprintf("previous position was close to j=%d, shifted point is not\n", j);
#endif
	    AUX(k,j) -= 1; /* point j loses a neighbour */
	  }
	}
      }
    }
#ifdef DEBUG
  Rprintf("end update ---- \n");
  for(j=0; j < npts; j++) {
    for(k=0; k < ndisc; k++)
      Rprintf("aux[%d,%d]=%d\t", k, j, AUX(k,j));
    Rprintf("\n");
  }
#endif
    return;
  }
  fexitc("Unrecognised transition type; bailing out.\n");
}

Cifns BadGeyCifns = { &badgeyinit, &badgeycif, &badgeyupd, FALSE};
