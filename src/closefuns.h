/*
  closefuns.h

  Function definitions to be #included in closepair.c
  several times with different values of macros.

  Macros used:

  CLOSEFUN   name of function for 'closepairs'

  CROSSFUN   name of function for 'crosspairs'

  COORDS     if defined, also return xi, yi, xj, yj, dx, dy, d

  THRESH     if defined, also return 1(d < s)

  ZCOORD     if defined, coordinates are 3-dimensional

  $Revision: 1.3 $ $Date: 2013/05/22 10:21:28 $

*/

#ifdef ZCOORD 
#define SPACEDIM 3
#else
#define SPACEDIM 2
#endif


SEXP CLOSEFUN(SEXP xx,
	      SEXP yy,
#ifdef ZCOORD
	      SEXP zz,
#endif
	      SEXP rr,
#ifdef THRESH
	      SEXP ss,
#endif
	      SEXP nguess) 
{
  double *x, *y;
  double xi, yi, rmax, r2max, dx, dy, d2;
#ifdef ZCOORD
  double *z;
  double zi, dz;
#endif
  int n, k, kmax, kmaxold, maxchunk, i, j, m;
  /* local storage */
  int *iout, *jout;
  /* R objects in return value */
  SEXP Out, iOut, jOut;
  /* external storage pointers */
  int *iOutP, *jOutP;

#ifdef COORDS
  double *xiout, *yiout, *xjout, *yjout, *dxout, *dyout, *dout;
  SEXP xiOut, yiOut, xjOut, yjOut, dxOut, dyOut, dOut;
  double *xiOutP, *yiOutP, *xjOutP, *yjOutP, *dxOutP, *dyOutP, *dOutP;
#ifdef ZCOORD
  double *ziout, *zjout, *dzout;
  SEXP ziOut, zjOut, dzOut;
  double *ziOutP, *zjOutP, *dzOutP;
#endif
#endif

#ifdef THRESH
  double s, s2;
  int *tout;
  SEXP tOut;
  int *tOutP;
#endif
  
  /* protect R objects from garbage collector */
  PROTECT(xx     = AS_NUMERIC(xx));
  PROTECT(yy     = AS_NUMERIC(yy));
#ifdef ZCOORD
  PROTECT(zz     = AS_NUMERIC(zz));
#endif
  PROTECT(rr     = AS_NUMERIC(rr));
  PROTECT(nguess = AS_INTEGER(nguess));
#ifdef THRESH
  PROTECT(ss     = AS_NUMERIC(ss));
#define NINPUTS (3+SPACEDIM)
#else
#define NINPUTS (2+SPACEDIM)
#endif

  /* Translate arguments from R to C */

  x = NUMERIC_POINTER(xx);
  y = NUMERIC_POINTER(yy);
#ifdef ZCOORD
  z = NUMERIC_POINTER(zz);
#endif

  n = LENGTH(xx);
  rmax = *(NUMERIC_POINTER(rr));
  kmax = *(INTEGER_POINTER(nguess));
  
  r2max = rmax * rmax;

#ifdef THRESH
  s = *(NUMERIC_POINTER(ss));
  s2 = s * s;
#endif

  k = 0;   /* k is the next available storage location 
              and also the current length of the list */ 

  if(n > 0 && kmax > 0) {
    /* allocate space */
    iout = (int *) R_alloc(kmax, sizeof(int));
    jout = (int *) R_alloc(kmax, sizeof(int));
#ifdef COORDS
    xiout =  (double *) R_alloc(kmax, sizeof(double));
    yiout =  (double *) R_alloc(kmax, sizeof(double));
    xjout =  (double *) R_alloc(kmax, sizeof(double));
    yjout =  (double *) R_alloc(kmax, sizeof(double));
    dxout =  (double *) R_alloc(kmax, sizeof(double));
    dyout =  (double *) R_alloc(kmax, sizeof(double));
#ifdef ZCOORD
    ziout =  (double *) R_alloc(kmax, sizeof(double));
    zjout =  (double *) R_alloc(kmax, sizeof(double));
    dzout =  (double *) R_alloc(kmax, sizeof(double));
#endif
    dout  =  (double *) R_alloc(kmax, sizeof(double));
#endif

#ifdef THRESH
    tout  =  (int *) R_alloc(kmax, sizeof(int));
#endif
    
    /* loop in chunks of 2^16 */

    i = 0; maxchunk = 0; 
    while(i < n) {

      R_CheckUserInterrupt();

      maxchunk += 65536; 
      if(maxchunk > n) maxchunk = n;

      for(; i < maxchunk; i++) {

	xi = x[i];
	yi = y[i];
#ifdef ZCOORD
	zi = z[i];
#endif
	if(i > 0) {
	  /* scan backward */
	  for(j = i - 1; j >= 0; j--) {
	    dx = x[j] - xi;
	    if(dx < -rmax) 
	      break;
	    dy = y[j] - yi;
	    d2 = dx * dx + dy * dy;
#ifdef ZCOORD
	    if(d2 <= r2max) {
	      dz = z[j] - zi;
	      d2 = d2 + dz * dz;
#endif
	      if(d2 <= r2max) {
		/* add this (i, j) pair to output */
		if(k >= kmax) {
		  /* overflow; allocate more space */
		  kmaxold = kmax;
		  kmax    = 2 * kmax;
		  iout  = intRealloc(iout,  kmaxold, kmax);
		  jout  = intRealloc(jout,  kmaxold, kmax);
#ifdef COORDS
		  xiout = dblRealloc(xiout, kmaxold, kmax); 
		  yiout = dblRealloc(yiout, kmaxold, kmax); 
		  xjout = dblRealloc(xjout, kmaxold, kmax); 
		  yjout = dblRealloc(yjout, kmaxold, kmax); 
		  dxout = dblRealloc(dxout, kmaxold, kmax); 
		  dyout = dblRealloc(dyout, kmaxold, kmax); 
#ifdef ZCOORD
		  ziout = dblRealloc(ziout, kmaxold, kmax); 
		  zjout = dblRealloc(zjout, kmaxold, kmax); 
		  dzout = dblRealloc(dzout, kmaxold, kmax); 
#endif
		  dout  = dblRealloc(dout,  kmaxold, kmax); 
#endif
#ifdef THRESH
		tout  = intRealloc(tout,  kmaxold, kmax);
#endif
	      }
	      jout[k] = j + 1; /* R indexing */
	      iout[k] = i + 1;
#ifdef COORDS
	      xiout[k] = xi;
	      yiout[k] = yi;
	      xjout[k] = x[j];
	      yjout[k] = y[j];
	      dxout[k] = dx;
	      dyout[k] = dy;
#ifdef ZCOORD
	      ziout[k] = zi;
	      zjout[k] = z[j];
	      dzout[k] = dz;
#endif
	      dout[k] = sqrt(d2);
#endif

#ifdef THRESH
	      tout[k] = (d2 <= s2) ? 1 : 0;
#endif
	      ++k;
	      }
#ifdef ZCOORD
	    }
#endif
	  }

	}
      
	if(i + 1 < n) {
	  /* scan forward */
	  for(j = i + 1; j < n; j++) {
	    dx = x[j] - xi;
	    if(dx > rmax) 
	      break;
	    dy = y[j] - yi;
	    d2 = dx * dx + dy * dy;
#ifdef ZCOORD
	    if(d2 <= r2max) {
	      dz = z[j] - zi;
	      d2 = d2 + dz * dz;
#endif
	      if(d2 <= r2max) {
		/* add this (i, j) pair to output */
		if(k >= kmax) {
		  /* overflow; allocate more space */
		  kmaxold = kmax;
		  kmax    = 2 * kmax;
		  iout  = intRealloc(iout,  kmaxold, kmax);
		  jout  = intRealloc(jout,  kmaxold, kmax);
#ifdef COORDS
		  xiout = dblRealloc(xiout, kmaxold, kmax); 
		  yiout = dblRealloc(yiout, kmaxold, kmax); 
		  xjout = dblRealloc(xjout, kmaxold, kmax); 
		  yjout = dblRealloc(yjout, kmaxold, kmax); 
		  dxout = dblRealloc(dxout, kmaxold, kmax); 
		  dyout = dblRealloc(dyout, kmaxold, kmax); 
#ifdef ZCOORD
		  ziout = dblRealloc(ziout, kmaxold, kmax); 
		  zjout = dblRealloc(zjout, kmaxold, kmax); 
		  dzout = dblRealloc(dzout, kmaxold, kmax); 
#endif
		  dout  = dblRealloc(dout,  kmaxold, kmax); 
#endif
#ifdef THRESH
		  tout  = intRealloc(tout,  kmaxold, kmax);
#endif
		}
		jout[k] = j + 1; /* R indexing */
		iout[k] = i + 1;
#ifdef COORDS
		xiout[k] = xi;
		yiout[k] = yi;
		xjout[k] = x[j];
		yjout[k] = y[j];
		dxout[k] = dx;
		dyout[k] = dy;
#ifdef ZCOORD
		ziout[k] = zi;
		zjout[k] = z[j];
		dzout[k] = dz;
#endif
		dout[k] = sqrt(d2);
#endif
#ifdef THRESH
		tout[k] = (d2 <= s2) ? 1 : 0;
#endif
		++k;
	      }
#ifdef ZCOORD
	    }
#endif
	  }
	}
	/* end of i loop */
      }
    }
  }

  /* return a list of vectors */
  PROTECT(iOut  = NEW_INTEGER(k));
  PROTECT(jOut  = NEW_INTEGER(k));
#ifdef COORDS
  PROTECT(xiOut = NEW_NUMERIC(k));
  PROTECT(yiOut = NEW_NUMERIC(k));
  PROTECT(xjOut = NEW_NUMERIC(k));
  PROTECT(yjOut = NEW_NUMERIC(k));
  PROTECT(dxOut = NEW_NUMERIC(k));
  PROTECT(dyOut = NEW_NUMERIC(k));
#ifdef ZCOORD
  PROTECT(ziOut = NEW_NUMERIC(k));
  PROTECT(zjOut = NEW_NUMERIC(k));
  PROTECT(dzOut = NEW_NUMERIC(k));
#endif
  PROTECT(dOut  = NEW_NUMERIC(k));
#endif
#ifdef THRESH
  PROTECT(tOut = NEW_INTEGER(k));
#endif
  if(k > 0) {
    iOutP  = INTEGER_POINTER(iOut);
    jOutP  = INTEGER_POINTER(jOut);
#ifdef COORDS
    xiOutP = NUMERIC_POINTER(xiOut);
    yiOutP = NUMERIC_POINTER(yiOut);
    xjOutP = NUMERIC_POINTER(xjOut);
    yjOutP = NUMERIC_POINTER(yjOut);
    dxOutP = NUMERIC_POINTER(dxOut);
    dyOutP = NUMERIC_POINTER(dyOut);
#ifdef ZCOORD
    ziOutP = NUMERIC_POINTER(ziOut);
    zjOutP = NUMERIC_POINTER(zjOut);
    dzOutP = NUMERIC_POINTER(dzOut);
#endif
    dOutP  = NUMERIC_POINTER(dOut);
#endif
#ifdef THRESH
    tOutP  = INTEGER_POINTER(tOut);
#endif
    for(m = 0; m < k; m++) {
      iOutP[m] = iout[m];
      jOutP[m] = jout[m];
#ifdef COORDS
      xiOutP[m] = xiout[m];
      yiOutP[m] = yiout[m];
      xjOutP[m] = xjout[m];
      yjOutP[m] = yjout[m];
      dxOutP[m] = dxout[m];
      dyOutP[m] = dyout[m];
#ifdef ZCOORD
      ziOutP[m] = ziout[m];
      zjOutP[m] = zjout[m];
      dzOutP[m] = dzout[m];
#endif
      dOutP[m]  = dout[m];
#endif
#ifdef THRESH
      tOutP[m]  = tout[m];
#endif
    }
  }

#define HEAD 2
#ifdef THRESH
#define MIDDLE 1
#else
#define MIDDLE 0
#endif
#ifdef COORDS
#define TAIL (3*SPACEDIM + 1)
#else 
#define TAIL 0
#endif

  PROTECT(Out   = NEW_LIST(HEAD+MIDDLE+TAIL));

  SET_VECTOR_ELT(Out, 0,  iOut);
  SET_VECTOR_ELT(Out, 1,  jOut);

#ifdef THRESH
  SET_VECTOR_ELT(Out, HEAD,  tOut);
#endif

#ifdef COORDS
#ifdef ZCOORD
  SET_VECTOR_ELT(Out, HEAD+MIDDLE,   xiOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+1, yiOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+2, ziOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+3, xjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+4, yjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+5, zjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+6, dxOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+7, dyOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+8, dzOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+9, dOut);
#else
  SET_VECTOR_ELT(Out, HEAD+MIDDLE,   xiOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+1, yiOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+2, xjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+3, yjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+4, dxOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+5, dyOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+6, dOut);
#endif
#endif

  UNPROTECT(NINPUTS+1+HEAD+MIDDLE+TAIL);    /* 1 is for 'Out' itself */

  return(Out);
}

#undef NINPUTS
#undef HEAD
#undef MIDDLE
#undef TAIL

SEXP CROSSFUN(SEXP xx1,
	      SEXP yy1,
#ifdef ZCOORD
	      SEXP zz1,
#endif
	      SEXP xx2,
	      SEXP yy2,
#ifdef ZCOORD
	      SEXP zz2,
#endif
	      SEXP rr,
#ifdef THRESH
	      SEXP ss,
#endif
	      SEXP nguess) 
{
  /* input vectors */
  double *x1, *y1, *x2, *y2;
#ifdef ZCOORD
  double *z1, *z2;
#endif
  /* lengths */
  int n1, n2, nout, noutmax, noutmaxold, maxchunk;
  /* distance parameter */
  double rmax, r2max;
  /* indices */
  int i, j, jleft, m;
  /* temporary values */
  double x1i, y1i, xleft, dx, dy, dx2, d2;
#ifdef ZCOORD
  double z1i, dz;
#endif
  /* local storage */
  int *iout, *jout;
  /* R objects in return value */
  SEXP Out, iOut, jOut;
  /* external storage pointers */
  int *iOutP, *jOutP;
#ifdef COORDS
  SEXP xiOut, yiOut, xjOut, yjOut, dxOut, dyOut, dOut;
  double *xiOutP, *yiOutP, *xjOutP, *yjOutP, *dxOutP, *dyOutP, *dOutP;
  double *xiout, *yiout, *xjout, *yjout, *dxout, *dyout, *dout;
#ifdef ZCOORD
  SEXP ziOut, zjOut, dzOut;
  double *ziOutP, *zjOutP, *dzOutP;
  double *ziout, *zjout, *dzout;
#endif
#endif  
#ifdef THRESH
  double s, s2;
  int *tout;
  SEXP tOut;
  int *tOutP;
#endif
  /* protect R objects from garbage collector */
  PROTECT(xx1     = AS_NUMERIC(xx1));
  PROTECT(yy1     = AS_NUMERIC(yy1));
  PROTECT(xx2     = AS_NUMERIC(xx2));
  PROTECT(yy2     = AS_NUMERIC(yy2));
#ifdef ZCOORD
  PROTECT(zz1     = AS_NUMERIC(zz1));
  PROTECT(zz2     = AS_NUMERIC(zz2));
#endif

  PROTECT(rr     = AS_NUMERIC(rr));
  PROTECT(nguess = AS_INTEGER(nguess));
#ifdef THRESH
  PROTECT(ss     = AS_NUMERIC(ss));
#define NINPUTS (2*SPACEDIM + 3)
#else
#define NINPUTS (2*SPACEDIM + 2)
#endif

  /* Translate arguments from R to C */

  x1 = NUMERIC_POINTER(xx1);
  y1 = NUMERIC_POINTER(yy1);
  x2 = NUMERIC_POINTER(xx2);
  y2 = NUMERIC_POINTER(yy2);
#ifdef ZCOORD
  z1 = NUMERIC_POINTER(zz1);
  z2 = NUMERIC_POINTER(zz2);
#endif

  n1 = LENGTH(xx1);
  n2 = LENGTH(xx2);
  rmax = *(NUMERIC_POINTER(rr));
  noutmax = *(INTEGER_POINTER(nguess));
  
  r2max = rmax * rmax;

#ifdef THRESH
  s = *(NUMERIC_POINTER(ss));
  s2 = s * s;
#endif

  nout = 0;   /* nout is the next available storage location 
		 and also the current length of the list */ 

  if(n1 > 0 && n2 > 0 && noutmax > 0) {
    /* allocate space */
    iout = (int *) R_alloc(noutmax, sizeof(int));
    jout = (int *) R_alloc(noutmax, sizeof(int));
#ifdef COORDS
    xiout =  (double *) R_alloc(noutmax, sizeof(double));
    yiout =  (double *) R_alloc(noutmax, sizeof(double));
    xjout =  (double *) R_alloc(noutmax, sizeof(double));
    yjout =  (double *) R_alloc(noutmax, sizeof(double));
    dxout =  (double *) R_alloc(noutmax, sizeof(double));
    dyout =  (double *) R_alloc(noutmax, sizeof(double));
#ifdef ZCOORD
    ziout =  (double *) R_alloc(noutmax, sizeof(double));
    zjout =  (double *) R_alloc(noutmax, sizeof(double));
    dzout =  (double *) R_alloc(noutmax, sizeof(double));
#endif
    dout  =  (double *) R_alloc(noutmax, sizeof(double));
#endif
#ifdef THRESH
    tout  =  (int *) R_alloc(noutmax, sizeof(int));
#endif
    
    jleft = 0;

    i = 0; maxchunk = 0;

    while(i < n1) {

      R_CheckUserInterrupt();

      maxchunk += 65536;
      if(maxchunk > n1) maxchunk = n1;

      for( ; i < maxchunk; i++) {

	x1i = x1[i];
	y1i = y1[i];
#ifdef ZCOORD
	z1i = z1[i];
#endif

	/* 
	   adjust starting point jleft
	*/
	xleft = x1i - rmax;
	while((x2[jleft] < xleft) && (jleft+1 < n2))
	  ++jleft;

	/* 
	   process from j = jleft until dx > rmax
	*/
	for(j=jleft; j < n2; j++) {

	  /* squared interpoint distance */
	  dx = x2[j] - x1i;
	  if(dx > rmax)
	    break;
	  dx2 = dx * dx;
	  dy = y2[j] - y1i;
	  d2 = dx2 + dy * dy;
#ifdef ZCOORD
	    if(d2 <= r2max) {
	      dz = z2[j] - z1i;
	      d2 = d2 + dz * dz;
#endif
	      if(d2 <= r2max) {
		/* add this (i, j) pair to output */
		if(nout >= noutmax) {
		  /* overflow; allocate more space */
		  noutmaxold = noutmax;
		  noutmax    = 2 * noutmax;
		  iout  = intRealloc(iout,  noutmaxold, noutmax);
		  jout  = intRealloc(jout,  noutmaxold, noutmax);
#ifdef COORDS
		  xiout = dblRealloc(xiout, noutmaxold, noutmax); 
		  yiout = dblRealloc(yiout, noutmaxold, noutmax); 
		  xjout = dblRealloc(xjout, noutmaxold, noutmax); 
		  yjout = dblRealloc(yjout, noutmaxold, noutmax); 
		  dxout = dblRealloc(dxout, noutmaxold, noutmax); 
		  dyout = dblRealloc(dyout, noutmaxold, noutmax); 
#ifdef ZCOORD
		  ziout = dblRealloc(ziout, noutmaxold, noutmax); 
		  zjout = dblRealloc(zjout, noutmaxold, noutmax); 
		  dzout = dblRealloc(dzout, noutmaxold, noutmax); 
#endif
		  dout  = dblRealloc(dout,  noutmaxold, noutmax); 
#endif
#ifdef THRESH
		  tout  = intRealloc(tout,  noutmaxold, noutmax);
#endif
		}
		iout[nout] = i + 1; /* R indexing */
		jout[nout] = j + 1;
#ifdef COORDS
		xiout[nout] = x1i;
		yiout[nout] = y1i;
		xjout[nout] = x2[j];
		yjout[nout] = y2[j];
		dxout[nout] = dx;
		dyout[nout] = dy;
#ifdef ZCOORD
		ziout[nout] = z1i;
		zjout[nout] = z2[j];
		dzout[nout] = dz;
#endif
		dout[nout] = sqrt(d2);
#endif
#ifdef THRESH
		tout[nout] = (d2 <= s2) ? 1 : 0;
#endif
		++nout;
	      }
#ifdef ZCOORD
	    }
#endif
	}
      }
    }
  }

  /* return a list of vectors */
  PROTECT(iOut  = NEW_INTEGER(nout));
  PROTECT(jOut  = NEW_INTEGER(nout));
#ifdef COORDS
  PROTECT(xiOut = NEW_NUMERIC(nout));
  PROTECT(yiOut = NEW_NUMERIC(nout));
  PROTECT(xjOut = NEW_NUMERIC(nout));
  PROTECT(yjOut = NEW_NUMERIC(nout));
  PROTECT(dxOut = NEW_NUMERIC(nout));
  PROTECT(dyOut = NEW_NUMERIC(nout));
#ifdef ZCOORD
  PROTECT(ziOut = NEW_NUMERIC(nout));
  PROTECT(zjOut = NEW_NUMERIC(nout));
  PROTECT(dzOut = NEW_NUMERIC(nout));
#endif
  PROTECT(dOut  = NEW_NUMERIC(nout));
#endif
#ifdef THRESH
  PROTECT(tOut = NEW_INTEGER(nout));
#endif
  if(nout > 0) {
    iOutP  = INTEGER_POINTER(iOut);
    jOutP  = INTEGER_POINTER(jOut);
#ifdef COORDS
    xiOutP = NUMERIC_POINTER(xiOut);
    yiOutP = NUMERIC_POINTER(yiOut);
    xjOutP = NUMERIC_POINTER(xjOut);
    yjOutP = NUMERIC_POINTER(yjOut);
    dxOutP = NUMERIC_POINTER(dxOut);
    dyOutP = NUMERIC_POINTER(dyOut);
#ifdef ZCOORD
    ziOutP = NUMERIC_POINTER(ziOut);
    zjOutP = NUMERIC_POINTER(zjOut);
    dzOutP = NUMERIC_POINTER(dzOut);
#endif
    dOutP  = NUMERIC_POINTER(dOut);
#endif
#ifdef THRESH
    tOutP  = INTEGER_POINTER(tOut);
#endif
    for(m = 0; m < nout; m++) {
      iOutP[m] = iout[m];
      jOutP[m] = jout[m];
#ifdef COORDS
      xiOutP[m] = xiout[m];
      yiOutP[m] = yiout[m];
      xjOutP[m] = xjout[m];
      yjOutP[m] = yjout[m];
      dxOutP[m] = dxout[m];
      dyOutP[m] = dyout[m];
#ifdef ZCOORD
      ziOutP[m] = ziout[m];
      zjOutP[m] = zjout[m];
      dzOutP[m] = dzout[m];
#endif
      dOutP[m]  = dout[m];
#endif
#ifdef THRESH
      tOutP[m]  = tout[m];
#endif
    }
  }
#define HEAD 2
#ifdef THRESH
#define MIDDLE 1
#else
#define MIDDLE 0
#endif
#ifdef COORDS
#define TAIL (3*SPACEDIM + 1)
#else 
#define TAIL 0
#endif

  PROTECT(Out   = NEW_LIST(HEAD+MIDDLE+TAIL));

  SET_VECTOR_ELT(Out, 0,  iOut);
  SET_VECTOR_ELT(Out, 1,  jOut);

#ifdef THRESH
  SET_VECTOR_ELT(Out, HEAD,  tOut);
#endif

#ifdef COORDS
#ifdef ZCOORD
  SET_VECTOR_ELT(Out, HEAD+MIDDLE,   xiOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+1, yiOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+2, ziOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+3, xjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+4, yjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+5, zjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+6, dxOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+7, dyOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+8, dzOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+9, dOut);
#else
  SET_VECTOR_ELT(Out, HEAD+MIDDLE,   xiOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+1, yiOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+2, xjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+3, yjOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+4, dxOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+5, dyOut);
  SET_VECTOR_ELT(Out, HEAD+MIDDLE+6, dOut);
#endif
#endif

  UNPROTECT(NINPUTS+1+HEAD+MIDDLE+TAIL);   /* 1 is for 'Out' itself */

  return(Out);
}

#undef NINPUTS
#undef HEAD
#undef MIDDLE
#undef TAIL



