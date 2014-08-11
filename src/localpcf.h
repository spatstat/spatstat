/*
  
  localpcf.h

  Source template for versions of local pair correlation

  Requires variable: WEIGHTED

  Assumes point patterns are sorted in increasing order of x coordinate

  $Revision: 1.5 $  $Date: 2012/03/27 04:50:04 $

*/

#ifdef WEIGHTED
#define FNAME locWpcfx
#else
#define FNAME locpcfx
#endif

void FNAME(nn1, x1, y1, id1, 
	   nn2, x2, y2, id2, 
#ifdef WEIGHTED
           w2,
#endif
	   nnr, rmaxi, 
	   del, pcf)
     /* inputs */
     int *nn1, *nn2, *nnr;
     double *x1, *y1, *x2, *y2;
     int *id1, *id2;
     double *rmaxi, *del;
#ifdef WEIGHTED
     double *w2;
#endif
     /* output */
     double *pcf;  /* matrix of column vectors of pcf's 
		      for each point of first pattern */
{
  int n1, n2, nr, i, j, k, jleft, kmin, kmax, id1i, maxchunk;
  double x1i, y1i, rmax, delta, xleft, dx, dy, dx2;
  double d2, d2max, dmax, d;
  double rstep, rvalue, frac, contrib, weight, coef;

  n1 = *nn1;
  n2 = *nn2;
  nr = *nnr;
  rmax = *rmaxi;
  delta = *del;

  dmax = rmax + delta; /* maximum relevant value of interpoint distance */
  d2max = dmax * dmax;
  rstep = rmax/(nr-1);
  coef  = 3.0 /(4.0 * delta);

  if(n1 == 0 || n2 == 0) 
    return;

  jleft = 0;

  OUTERCHUNKLOOP(i, n1, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, n1, maxchunk, 8196) {
      x1i = x1[i];
      y1i = y1[i];
      id1i = id1[i];

      /* 
	 adjust starting point

      */
      xleft = x1i - dmax;
      while((x2[jleft] < xleft) && (jleft+1 < n2))
	++jleft;

      /* 
	 process from jleft until |dx| > dmax
      */
      for(j=jleft; j < n2; j++) {
	dx = x2[j] - x1i;
	dx2 = dx * dx;
	if(dx2 > d2max) 
	  break;
	dy = y2[j] - y1i;
	d2 = dx2 + dy * dy;
	if(d2 <= d2max && id2[j] != id1i) {
	  d = sqrt(d2);
	  kmin = (int) floor((d-delta)/rstep);
	  kmax = (int) ceil((d+delta)/rstep);
	  if(kmin <= nr-1 && kmax >= 0) {
	    /* nonempty intersection with range of r values */
	    /* compute intersection */
	    if(kmin < 0) kmin = 0; 
	    if(kmax >= nr) kmax = nr-1;
	    /* */
	    weight = coef/d;
#ifdef WEIGHTED
	    weight = weight * w2[j];
#endif
	    for(k = kmin; k <= kmax; k++) {
	      rvalue = k * rstep;
	      frac = (d - rvalue)/delta;
	      /* Epanechnikov kernel with halfwidth delta */
	      contrib = (1 - frac * frac);
	      if(contrib > 0) 
		pcf[k + nr * i] += contrib * weight;
	    }
	  }
	}
      }
    }
  }
}

#undef FNAME
