#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "geom3.h"
#include "functable.h"
#include "chunkloop.h"
#include "constants.h"

/*
	$Revision: 1.7 $	$Date: 2012/03/27 05:01:41 $

	pair correlation function of 3D point pattern
	(Epanechnikov kernel) 

	pcf3trans	  	translation correction

	pcf3isot		isotropic correction

*/

#define FOURPI (2.0 * M_2PI)


void
pcf3trans(p, n, b, pcf, delta)
     Point *p;
     int n;
     Box *b;
     Ftable *pcf;
     double delta;
{
  register int i, j, l, lmin, lmax, maxchunk;
  register double dx, dy, dz, dist;
  register double  vx, vy, vz, tval;
  Point *ip, *jp;
  double dt, vol, lambda, denom;
  double coef, twocoef, frac, invweight, kernel;

  double sphesfrac(), sphevol();

  /* compute denominator & initialise numerator*/
  vol = (b->x1 - b->x0) * (b->y1 - b->y0) * (b->z1 - b->z0);
  lambda = ((double) n )/ vol;
  denom = lambda * lambda;

  for(l = 0; l < pcf->n; l++) {
    (pcf->denom)[l] = denom;
    (pcf->num)[l]   = 0.0;
  }

  /* spacing of argument in result vector */
  dt = (pcf->t1 - pcf->t0)/(pcf->n - 1);

  /* compute numerator */
  OUTERCHUNKLOOP(i, n, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, n, maxchunk, 8196) {
      ip = p + i;
      for(j = i + 1; j < n; j++) {
	/* compute pairwise distance */
	jp = p + j;
	dx = jp->x - ip->x;
	dy = jp->y - ip->y;
	dz = jp->z - ip->z;
	dist = sqrt(dx * dx + dy * dy + dz * dz);
	lmin = ceil( ((dist - delta) - pcf->t0) / dt );
	lmax = floor( ((dist + delta) - pcf->t0) / dt );
	if(lmax >= 0 && lmin < pcf->n) {
	  /* kernel centred at 'dist' has nonempty intersection 
	     with specified range of t values */
	  /* compute intersection */
	  if(lmin < 0)
	    lmin = 0;
	  if(lmax >= pcf->n)
	    lmax = pcf->n - 1;
	  /* compute (inverse) edge correction weight */
	  vx = b->x1 - b->x0 - (dx > 0 ? dx : -dx);
	  vy = b->y1 - b->y0 - (dy > 0 ? dy : -dy);
	  vz = b->z1 - b->z0 - (dz > 0 ? dz : -dz);
	  invweight = vx * vy * vz * FOURPI * dist * dist;
	  if(invweight > 0.0) {
	    for(l = lmin; l < pcf->n; l++) {
	      tval = pcf->t0 + l * dt;
	      /* unnormalised Epanechnikov kernel with halfwidth delta */
	      frac = (dist - tval)/delta;
	      kernel = (1 - frac * frac);
	      if(kernel > 0) 	    
		(pcf->num)[l] += kernel / invweight;
	    }
	  }
	}
      }
    }
  }
  
  /* constant factor in kernel */
  coef = 3.0/(4.0 * delta);
  /* multiplied by 2 because we only visited i < j pairs */
  twocoef = 2.0 * coef; 

  /* normalise kernel and compute ratio estimate */
  for(l = 0; l < pcf->n; l++) {
    (pcf->num)[l] *= twocoef;
    (pcf->f)[l] = ((pcf->denom)[l] > 0.0) ?
      (pcf->num)[l] / (pcf->denom)[l] : 0.0;
  }
}


void
pcf3isot(p, n, b, pcf, delta)
     Point *p;
     int n;
     Box *b;
     Ftable *pcf;
     double delta;
{
  register int i, j, l, lmin, lmax, maxchunk;
  register double dx, dy, dz, dist;
  Point *ip, *jp;
  double dt, vol, denom, mass, tval;
  double coef, frac, kernel;

  double sphesfrac(), sphevol();
  Point vertex;
  Box   half;

  /* compute denominator & initialise numerator*/
  vol = (b->x1 - b->x0) * (b->y1 - b->y0) * (b->z1 - b->z0);
  denom = ((double) (n * n))/vol;

  for(l = 0; l < pcf->n; l++) {
    (pcf->denom)[l] = denom;
    (pcf->num)[l]   = 0.0;
  }

  /* spacing of argument in result vector */
  dt = (pcf->t1 - pcf->t0)/(pcf->n - 1);

  /* set up for volume correction */

  vertex.x = b->x0;
  vertex.y = b->y0;
  vertex.z = b->z0;
  half.x1  = b->x1;
  half.y1  = b->y1;
  half.z1  = b->z1;
  half.x0  = (b->x0 + b->x1)/2.0;
  half.y0  = (b->y0 + b->y1)/2.0;
  half.z0  = (b->z0 + b->z1)/2.0;

	/* compute numerator */
  OUTERCHUNKLOOP(i, n, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, n, maxchunk, 8196) {
      ip = p + i;
      for(j = i + 1; j < n; j++) {
	jp = p + j;
	dx = jp->x - ip->x;
	dy = jp->y - ip->y;
	dz = jp->z - ip->z;
	dist = sqrt(dx * dx + dy * dy + dz * dz);
	lmin = ceil( ((dist - delta) - pcf->t0) / dt );
	lmax = floor( ((dist + delta) - pcf->t0) / dt );
	if(lmax >= 0 && lmin < pcf->n) {
	  /* kernel centred at 'dist' has nonempty intersection 
	     with specified range of t values */
	  /* compute intersection */
	  if(lmin < 0)
	    lmin = 0;
	  if(lmax >= pcf->n)
	    lmax = pcf->n - 1;
	  /* compute edge correction weight */
	  mass = (1.0 / sphesfrac(ip, b, dist)) 
	    + (1.0 / sphesfrac(jp, b, dist)); 
	  mass *= 
	    1.0 - 8.0 * sphevol(&vertex, &half, dist) / vol;
	  if(mass > 0.0) {
	    mass /= FOURPI * dist * dist;
	    for(l = lmin; l < pcf->n; l++) {
	      tval = pcf->t0 + l * dt;
	      /* unnormalised Epanechnikov kernel with halfwidth delta */
	      frac = (dist - tval)/delta;
	      kernel = (1 - frac * frac);
	      if(kernel > 0) 	    
		(pcf->num)[l] += kernel * mass;
	    }
	  }
	}
      }
    }
  }

  /* constant factor in kernel */
  coef = 3.0/(4.0 * delta);

  /* normalise kernel and compute ratio estimate */
  for(l = 0; l < pcf->n; l++) {
    (pcf->num)[l] *= coef;
    (pcf->f)[l] = ((pcf->denom)[l] > 0.0)?
      (pcf->num)[l] / (pcf->denom)[l]
      : 0.0;
  }
}
