#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Triplets process */

/* Format for storage of parameters and precomputed/auxiliary data */

typedef struct Triplets {
  double gamma;
  double r;
  double loggamma;
  double r2;
  double *period;
  int hard;
  int per;
  int *neighbour;    /* scratch list of neighbours of current point */
  int Nmax;          /* length of scratch space allocated */
} Triplets;

/* initialiser function */

Cdata *tripletsinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* create storage for model parameters */
  Triplets *triplets;
  triplets = (Triplets *) R_alloc(1, sizeof(Triplets)); 
  /* create scratch space */
  triplets->Nmax  = 1024;
  triplets->neighbour = (int *) R_alloc(1024, sizeof(int));
  /* Interpret model parameters*/
  triplets->gamma  = model.ipar[0];
  triplets->r      = model.ipar[1]; /* No longer passed as r^2 */
  triplets->r2     = triplets->r * triplets->r; 
  triplets->period = model.period;
#ifdef MHDEBUG
  Rprintf("Initialising Triplets gamma=%lf, r=%lf\n", 
	  triplets->gamma, triplets->r);
#endif
  /* is the model numerically equivalent to hard core ? */
  triplets->hard   = (triplets->gamma < DOUBLE_EPS);
  triplets->loggamma = (triplets->hard) ? 0 : log(triplets->gamma);
  /* periodic boundary conditions? */
  triplets->per    = (model.period[0] > 0.0);
  return((Cdata *) triplets);
}

/* conditional intensity evaluator */

double tripletscif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, kount, ix, ixp1, j, k, nj, nk, N, Nmax, Nmore, N1;
  int *neighbour;
  double *x, *y;
  double u, v;
  double r2, d2,  cifval;
  Triplets *triplets;

  triplets = (Triplets *) cdata;

  r2     = triplets->r2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;

  if(npts == 0) 
    return((double) 1.0);

  neighbour = triplets->neighbour;
  Nmax      = triplets->Nmax;
  N         = 0;

  /* compile list of neighbours */

  for(j=0; j < npts; j++) {
    if(j != ix) {
      d2 = dist2either(u,v,x[j],y[j],triplets->period);
      if(d2 < r2) {
	/* add j to list of neighbours of current point */
	if(N >= Nmax) {
	  /* storage space overflow: reallocate */
	  Nmore = 2 * Nmax;
	  triplets->neighbour = neighbour = 
	    (int *) S_realloc((char *) triplets->neighbour,
			      Nmore, Nmax, sizeof(int));
	  triplets->Nmax = Nmax = Nmore;
	}
	neighbour[N] = j;
	N++;
      }
    }
  }

  /* count r-close (ordered) pairs of neighbours */
  kount = 0;

  if(N > 1) {
    N1 = N - 1;
    for(j = 0; j < N1; j++) {
      nj = neighbour[j];
      for(k = j+1; k < N; k++) {
	nk = neighbour[k];
	if(nj != nk) {
	  d2 = dist2either(x[nj],y[nj],x[nk],y[nk],triplets->period);
	  if(d2 < r2) kount++;
	}
      }
    }
  }
  
if(triplets->hard) {
    if(kount > 0) cifval = 0.0;
    else cifval = 1.0;
  } else cifval = exp((triplets->loggamma) * kount);

#ifdef MHDEBUG
 Rprintf("triplet count=%d cif=%lf\n", kount, cifval);
#endif

  return cifval;
}

Cifns TripletsCifns = { &tripletsinit, &tripletscif, (updafunptr) NULL, NO};
