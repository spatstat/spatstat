/*
  closefuns.h

  Function definitions to be #included in closepair.c
  several times with different values of macros.

  Macros used:

  CLOSEFUN   name of function for 'closepairs'

  CROSSFUN   name of function for 'crosspairs'

  EVERYTHING  if not defined, return i, j 
              if defined,     return i, j, xi, yi, xj, yj, dx, dy, d

  $Revision: 1.1 $ $Date: 2013/02/22 01:05:42 $

*/


SEXP CLOSEFUN(SEXP xx,
		 SEXP yy,
		 SEXP rr,
		 SEXP nguess) 
{
  double *x, *y;
  double xi, yi, rmax, r2max, dx, dy, dx2, d2;
  int n, k, kmax, kmaxold, maxchunk, i, j, m;
  /* local storage */
  int *iout, *jout;
  /* R objects in return value */
  SEXP Out, iOut, jOut;
  /* external storage pointers */
  int *iOutP, *jOutP;
#ifdef EVERYTHING
  double *xiout, *yiout, *xjout, *yjout, *dxout, *dyout, *dout;
  SEXP xiOut, yiOut, xjOut, yjOut, dxOut, dyOut, dOut;
  double *xiOutP, *yiOutP, *xjOutP, *yjOutP, *dxOutP, *dyOutP, *dOutP;
#endif
  
  /* protect R objects from garbage collector */
  PROTECT(xx     = AS_NUMERIC(xx));
  PROTECT(yy     = AS_NUMERIC(yy));
  PROTECT(rr     = AS_NUMERIC(rr));
  PROTECT(nguess = AS_INTEGER(nguess));
  /* That's 4 objects */

  /* Translate arguments from R to C */

  x = NUMERIC_POINTER(xx);
  y = NUMERIC_POINTER(yy);
  n = LENGTH(xx);
  rmax = *(NUMERIC_POINTER(rr));
  kmax = *(INTEGER_POINTER(nguess));
  
  r2max = rmax * rmax;

  k = 0;   /* k is the next available storage location 
              and also the current length of the list */ 

  if(n > 0 && kmax > 0) {
    /* allocate space */
    iout = (int *) R_alloc(kmax, sizeof(int));
    jout = (int *) R_alloc(kmax, sizeof(int));
#ifdef EVERYTHING
    xiout =  (double *) R_alloc(kmax, sizeof(double));
    yiout =  (double *) R_alloc(kmax, sizeof(double));
    xjout =  (double *) R_alloc(kmax, sizeof(double));
    yjout =  (double *) R_alloc(kmax, sizeof(double));
    dxout =  (double *) R_alloc(kmax, sizeof(double));
    dyout =  (double *) R_alloc(kmax, sizeof(double));
    dout  =  (double *) R_alloc(kmax, sizeof(double));
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

	if(i > 0) {
	  /* scan backward */
	  for(j = i - 1; j >= 0; j--) {
	    dx = x[j] - xi;
	    dx2 = dx * dx;
	    if(dx2 > r2max)
	      break;
	    dy = y[j] - yi;
	    d2 = dx2 + dy * dy;
	    if(d2 <= r2max) {
	      /* add this (i, j) pair to output */
	      if(k >= kmax) {
		/* overflow; allocate more space */
		kmaxold = kmax;
		kmax    = 2 * kmax;
		iout  = intRealloc(iout,  kmaxold, kmax);
		jout  = intRealloc(jout,  kmaxold, kmax);
#ifdef EVERYTHING
		xiout = dblRealloc(xiout, kmaxold, kmax); 
		yiout = dblRealloc(yiout, kmaxold, kmax); 
		xjout = dblRealloc(xjout, kmaxold, kmax); 
		yjout = dblRealloc(yjout, kmaxold, kmax); 
		dxout = dblRealloc(dxout, kmaxold, kmax); 
		dyout = dblRealloc(dyout, kmaxold, kmax); 
		dout  = dblRealloc(dout,  kmaxold, kmax); 
#endif
	      }
	      jout[k] = j + 1; /* R indexing */
	      iout[k] = i + 1;
#ifdef EVERYTHING
	      xiout[k] = xi;
	      yiout[k] = yi;
	      xjout[k] = x[j];
	      yjout[k] = y[j];
	      dxout[k] = dx;
	      dyout[k] = dy;
	      dout[k] = sqrt(d2);
#endif
	      ++k;
	    }
	  }
	}
      
	if(i + 1 < n) {
	  /* scan forward */
	  for(j = i + 1; j < n; j++) {
	    dx = x[j] - xi;
	    dx2 = dx * dx;
	    if(dx2 > r2max)
	      break;
	    dy = y[j] - yi;
	    d2 = dx2 + dy * dy;
	    if(d2 <= r2max) {
	      /* add this (i, j) pair to output */
	      if(k >= kmax) {
		/* overflow; allocate more space */
		kmaxold = kmax;
		kmax    = 2 * kmax;
		iout  = intRealloc(iout,  kmaxold, kmax);
		jout  = intRealloc(jout,  kmaxold, kmax);
#ifdef EVERYTHING
		xiout = dblRealloc(xiout, kmaxold, kmax); 
		yiout = dblRealloc(yiout, kmaxold, kmax); 
		xjout = dblRealloc(xjout, kmaxold, kmax); 
		yjout = dblRealloc(yjout, kmaxold, kmax); 
		dxout = dblRealloc(dxout, kmaxold, kmax); 
		dyout = dblRealloc(dyout, kmaxold, kmax); 
		dout  = dblRealloc(dout,  kmaxold, kmax); 
#endif
	      }
	      jout[k] = j + 1; /* R indexing */
	      iout[k] = i + 1;
#ifdef EVERYTHING
	      xiout[k] = xi;
	      yiout[k] = yi;
	      xjout[k] = x[j];
	      yjout[k] = y[j];
	      dxout[k] = dx;
	      dyout[k] = dy;
	      dout[k] = sqrt(d2);
#endif
	      ++k;
	    }
	  }
	}
	/* end of i loop */
      }
    }
  }

  /* return a list of vectors */
  PROTECT(iOut  = NEW_INTEGER(k));
  PROTECT(jOut  = NEW_INTEGER(k));
#ifdef EVERYTHING
  PROTECT(xiOut = NEW_NUMERIC(k));
  PROTECT(yiOut = NEW_NUMERIC(k));
  PROTECT(xjOut = NEW_NUMERIC(k));
  PROTECT(yjOut = NEW_NUMERIC(k));
  PROTECT(dxOut = NEW_NUMERIC(k));
  PROTECT(dyOut = NEW_NUMERIC(k));
  PROTECT(dOut  = NEW_NUMERIC(k));
#endif
  if(k > 0) {
    iOutP  = INTEGER_POINTER(iOut);
    jOutP  = INTEGER_POINTER(jOut);
#ifdef EVERYTHING
    xiOutP = NUMERIC_POINTER(xiOut);
    yiOutP = NUMERIC_POINTER(yiOut);
    xjOutP = NUMERIC_POINTER(xjOut);
    yjOutP = NUMERIC_POINTER(yjOut);
    dxOutP = NUMERIC_POINTER(dxOut);
    dyOutP = NUMERIC_POINTER(dyOut);
    dOutP  = NUMERIC_POINTER(dOut);
#endif
    for(m = 0; m < k; m++) {
      iOutP[m] = iout[m];
      jOutP[m] = jout[m];
#ifdef EVERYTHING
      xiOutP[m] = xiout[m];
      yiOutP[m] = yiout[m];
      xjOutP[m] = xjout[m];
      yjOutP[m] = yjout[m];
      dxOutP[m] = dxout[m];
      dyOutP[m] = dyout[m];
      dOutP[m]  = dout[m];
#endif
    }
  }
#ifdef EVERYTHING
  PROTECT(Out   = NEW_LIST(9));
#else 
  PROTECT(Out   = NEW_LIST(2));
#endif
  SET_VECTOR_ELT(Out, 0,  iOut);
  SET_VECTOR_ELT(Out, 1,  jOut);
#ifdef EVERYTHING
  SET_VECTOR_ELT(Out, 2, xiOut);
  SET_VECTOR_ELT(Out, 3, yiOut);
  SET_VECTOR_ELT(Out, 4, xjOut);
  SET_VECTOR_ELT(Out, 5, yjOut);
  SET_VECTOR_ELT(Out, 6, dxOut);
  SET_VECTOR_ELT(Out, 7, dyOut);
  SET_VECTOR_ELT(Out, 8, dOut);
#endif
#ifdef EVERYTHING
  UNPROTECT(14); /* 4 inputs and 10 outputs (Out and its 9 components) */
#else 
  UNPROTECT(7);  /* 4 inputs and 3 outputs (Out and its 2 components) */
#endif
  return(Out);
}



SEXP CROSSFUN(SEXP xx1,
	      SEXP yy1,
	      SEXP xx2,
	      SEXP yy2,
	      SEXP rr,
	      SEXP nguess) 
{
  /* input vectors */
  double *x1, *y1, *x2, *y2;
  /* lengths */
  int n1, n2, nout, noutmax, noutmaxold, maxchunk;
  /* distance parameter */
  double rmax, r2max;
  /* indices */
  int i, j, jleft, m;
  /* temporary values */
  double x1i, y1i, xleft, dx, dy, dx2, d2;
  /* local storage */
  int *iout, *jout;
  /* R objects in return value */
  SEXP Out, iOut, jOut;
  /* external storage pointers */
  int *iOutP, *jOutP;
#ifdef EVERYTHING
  SEXP xiOut, yiOut, xjOut, yjOut, dxOut, dyOut, dOut;
  double *xiOutP, *yiOutP, *xjOutP, *yjOutP, *dxOutP, *dyOutP, *dOutP;
  double *xiout, *yiout, *xjout, *yjout, *dxout, *dyout, *dout;
#endif  
  /* protect R objects from garbage collector */
  PROTECT(xx1     = AS_NUMERIC(xx1));
  PROTECT(yy1     = AS_NUMERIC(yy1));
  PROTECT(xx2     = AS_NUMERIC(xx2));
  PROTECT(yy2     = AS_NUMERIC(yy2));
  PROTECT(rr     = AS_NUMERIC(rr));
  PROTECT(nguess = AS_INTEGER(nguess));
  /* That's 6 objects */

  /* Translate arguments from R to C */

  x1 = NUMERIC_POINTER(xx1);
  y1 = NUMERIC_POINTER(yy1);
  x2 = NUMERIC_POINTER(xx2);
  y2 = NUMERIC_POINTER(yy2);
  n1 = LENGTH(xx1);
  n2 = LENGTH(xx2);
  rmax = *(NUMERIC_POINTER(rr));
  noutmax = *(INTEGER_POINTER(nguess));
  
  r2max = rmax * rmax;

  nout = 0;   /* nout is the next available storage location 
		 and also the current length of the list */ 

  if(n1 > 0 && n2 > 0 && noutmax > 0) {
    /* allocate space */
    iout = (int *) R_alloc(noutmax, sizeof(int));
    jout = (int *) R_alloc(noutmax, sizeof(int));
#ifdef EVERYTHING
    xiout =  (double *) R_alloc(noutmax, sizeof(double));
    yiout =  (double *) R_alloc(noutmax, sizeof(double));
    xjout =  (double *) R_alloc(noutmax, sizeof(double));
    yjout =  (double *) R_alloc(noutmax, sizeof(double));
    dxout =  (double *) R_alloc(noutmax, sizeof(double));
    dyout =  (double *) R_alloc(noutmax, sizeof(double));
    dout  =  (double *) R_alloc(noutmax, sizeof(double));
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
	  dx2 = dx * dx;
	  if(dx2 > r2max)
	    break;
	  dy = y2[j] - y1i;
	  d2 = dx2 + dy * dy;
	  if(d2 <= r2max) {
	    /* add this (i, j) pair to output */
	    if(nout >= noutmax) {
	      /* overflow; allocate more space */
	      noutmaxold = noutmax;
	      noutmax    = 2 * noutmax;
	      iout  = intRealloc(iout,  noutmaxold, noutmax);
	      jout  = intRealloc(jout,  noutmaxold, noutmax);
#ifdef EVERYTHING
	      xiout = dblRealloc(xiout, noutmaxold, noutmax); 
	      yiout = dblRealloc(yiout, noutmaxold, noutmax); 
	      xjout = dblRealloc(xjout, noutmaxold, noutmax); 
	      yjout = dblRealloc(yjout, noutmaxold, noutmax); 
	      dxout = dblRealloc(dxout, noutmaxold, noutmax); 
	      dyout = dblRealloc(dyout, noutmaxold, noutmax); 
	      dout  = dblRealloc(dout,  noutmaxold, noutmax); 
#endif
	    }
	    iout[nout] = i + 1; /* R indexing */
	    jout[nout] = j + 1;
#ifdef EVERYTHING
	    xiout[nout] = x1i;
	    yiout[nout] = y1i;
	    xjout[nout] = x2[j];
	    yjout[nout] = y2[j];
	    dxout[nout] = dx;
	    dyout[nout] = dy;
	    dout[nout] = sqrt(d2);
#endif
	    ++nout;
	  }
	}
      }
    }
  }

  /* return a list of vectors */
  PROTECT(iOut  = NEW_INTEGER(nout));
  PROTECT(jOut  = NEW_INTEGER(nout));
#ifdef EVERYTHING
  PROTECT(xiOut = NEW_NUMERIC(nout));
  PROTECT(yiOut = NEW_NUMERIC(nout));
  PROTECT(xjOut = NEW_NUMERIC(nout));
  PROTECT(yjOut = NEW_NUMERIC(nout));
  PROTECT(dxOut = NEW_NUMERIC(nout));
  PROTECT(dyOut = NEW_NUMERIC(nout));
  PROTECT(dOut  = NEW_NUMERIC(nout));
#endif
  if(nout > 0) {
    iOutP  = INTEGER_POINTER(iOut);
    jOutP  = INTEGER_POINTER(jOut);
#ifdef EVERYTHING
    xiOutP = NUMERIC_POINTER(xiOut);
    yiOutP = NUMERIC_POINTER(yiOut);
    xjOutP = NUMERIC_POINTER(xjOut);
    yjOutP = NUMERIC_POINTER(yjOut);
    dxOutP = NUMERIC_POINTER(dxOut);
    dyOutP = NUMERIC_POINTER(dyOut);
    dOutP  = NUMERIC_POINTER(dOut);
#endif
    for(m = 0; m < nout; m++) {
      iOutP[m] = iout[m];
      jOutP[m] = jout[m];
#ifdef EVERYTHING
      xiOutP[m] = xiout[m];
      yiOutP[m] = yiout[m];
      xjOutP[m] = xjout[m];
      yjOutP[m] = yjout[m];
      dxOutP[m] = dxout[m];
      dyOutP[m] = dyout[m];
      dOutP[m]  = dout[m];
#endif
    }
  }
#ifdef EVERYTHING
  PROTECT(Out   = NEW_LIST(9));
#else
  PROTECT(Out   = NEW_LIST(2));
#endif
  SET_VECTOR_ELT(Out, 0,  iOut);
  SET_VECTOR_ELT(Out, 1,  jOut);
#ifdef EVERYTHING
  SET_VECTOR_ELT(Out, 2, xiOut);
  SET_VECTOR_ELT(Out, 3, yiOut);
  SET_VECTOR_ELT(Out, 4, xjOut);
  SET_VECTOR_ELT(Out, 5, yjOut);
  SET_VECTOR_ELT(Out, 6, dxOut);
  SET_VECTOR_ELT(Out, 7, dyOut);
  SET_VECTOR_ELT(Out, 8, dOut);
#endif
#ifdef EVERYTHING
  UNPROTECT(16); /* 6 inputs and 10 outputs (Out and its 9 components) */
#else
  UNPROTECT(9); /* 6 inputs and 3 outputs (Out and its 2 components) */
#endif
  return(Out);
}



