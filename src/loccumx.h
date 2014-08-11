/*

  loccumx.h

  C template for loccum.c

  grid-to-data functions

  $Revision: 1.5 $ $Date: 2012/11/10 06:13:52 $

  macros: 

  FNAME    function name
  NULVAL   initial value (empty sum = 0, empty product = 1)
  INC(A,B) increment operation A += B or A *= B

*/

void FNAME(ntest, xtest, ytest, 
	   ndata, xdata, ydata, vdata,
	   nr, rmax, 
	   ans)
     /* inputs */
     int *ntest, *ndata, *nr;
     double *xtest, *ytest, *xdata, *ydata, *vdata;
     double *rmax;
     /* output */
     double *ans;  /* matrix of column vectors of functions 
		      for each point of first pattern */
{
  int Ntest, Ndata, Nr, Nans;
  double Rmax;

  int i, j, k, jleft, kmin, maxchunk, columnstart;
  double Rmax2, rstep, xtesti, ytesti, xleft;
  double dx, dy, dx2, d2, d, contrib;

  Ntest = *ntest;
  Ndata = *ndata;
  Nr    = *nr;
  Rmax  = *rmax;

  if(Ntest == 0)
    return;

  Nans  = Nr * Ntest;

  /* initialise products to 1 */
  OUTERCHUNKLOOP(k, Nans, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(k, Nans, maxchunk, 8196) {
      ans[k] = NULVAL;
    }
  }
   
  if(Ndata == 0) 
    return;

  rstep = Rmax/(Nr-1);
  Rmax2 = Rmax * Rmax;

  jleft = 0;

  OUTERCHUNKLOOP(i, Ntest, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Ntest, maxchunk, 8196) {
      xtesti = xtest[i];
      ytesti = ytest[i];
      columnstart = Nr * i; /* start position for f_i(.) in 'ans' */
      /* 
	 adjust starting point

      */
      xleft = xtesti - Rmax;
      while((xdata[jleft] < xleft) && (jleft+1 < Ndata))
	++jleft;

      /* 
	 process from jleft until |dx| > Rmax
      */
      for(j=jleft; j < Ndata; j++) {
	dx = xdata[j] - xtesti;
	dx2 = dx * dx;
	if(dx2 > Rmax2) 
	  break;
	dy = ydata[j] - ytesti;
	d2 = dx2 + dy * dy;
	if(d2 <= Rmax2) {
	  d = sqrt(d2);
	  kmin = (int) ceil(d/rstep);
	  if(kmin < Nr) {
	    contrib = vdata[j];
	    for(k = kmin; k < Nr; k++) 
	      INC(ans[columnstart + k] , contrib);
	  }
	}
      }
    }
  }
}

