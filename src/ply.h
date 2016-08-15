/*
    ply.h

    Template for functions in ply.c
    This file is #included several times

    Macros used:
    FNAME     Name of C routine
    NDIM      Number of dimensions of result (2 or 3)

    Adrian Baddeley and Tilman Davies

    $Revision: 1.1 $  $Date: 2016/08/15 02:29:15 $

*/


void FNAME(nin, 
	   xin,  
	   iin,  
	   jin,
#if (NDIM > 2)
	   kin,  
#endif
	   nout,
	   xout, 
	   iout, 
	   jout
#if (NDIM > 2)
	   , kout
#endif
) 
     int *nin, *nout;
     double *xin, *xout;
     int *iin, *jin, *iout, *jout;
#if (NDIM > 2)
     int *kin, *kout;
#endif
{
  int Nin, l, m, icur, jcur;
#if (NDIM > 2)
  int kcur;
#endif
  double xsum;
  Nin = *nin;
  if(Nin == 0) {
    *nout = 0;
    return;
  }
  /* initialise first cell using first entry */
  m = 0;
  iout[0] = icur = iin[0];
  jout[0] = jcur = jin[0];
#if (NDIM > 2)
  kout[0] = kcur = kin[0];
#endif
  xout[0] = xsum = xin[0];
  /* process subsequent entries */
  if(Nin > 1) {
    for(l = 1; l < Nin; l++) {
      if(iin[l] == icur && jin[l] == jcur 
#if (NDIM > 2)
	 && kin[l] == kcur
#endif
	 ) {
	/* increment current sum */
	xsum += xin[l];
      } else {
	/* write cell result */
	xout[m] = xsum;
	/* initialise next cell */
	++m;
	iout[m] = icur = iin[l];
	jout[m] = jcur = jin[l];
#if (NDIM > 2)
	kout[m] = kcur = kin[l];
#endif
	xsum = xin[l];
      }
      /* write last cell */
      xout[m] = xsum;
    }
  }
  *nout = m + 1;
}
