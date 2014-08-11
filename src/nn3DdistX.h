/*

  nn3DdistX.h

  Code template for nearest-neighbour algorithms for 3D point patterns

  Input is two point patterns - supports 'nncross'

  This code is #included multiple times in nn3Ddist.c
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
	EXCLUDE   #defined if the two patterns may include common points
	          (which are not to be counted as neighbours)

  Either or both DIST and WHICH may be defined.

  THE FOLLOWING CODE ASSUMES THAT BOTH POINT PATTERNS ARE SORTED
  IN ASCENDING ORDER OF THE z COORDINATE

  If EXCLUDE is #defined, 
   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].

  $Revision: 1.3 $ $Date: 2013/06/28 10:38:34 $

*/

void FNAME(n1, x1, y1, z1, id1, 
	   n2, x2, y2, z2, id2,
	   nnd, nnwhich, huge)
/* inputs */
     int *n1, *n2, *id1, *id2;
     double *x1, *y1, *z1, *x2, *y2, *z2, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int npoints1, npoints2, i, j, jwhich, lastjwhich, id1i, maxchunk;
  double dmin, d2, d2min, x1i, y1i, z1i, dx, dy, dz, dz2, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  lastjwhich = 0;

  for(i = 0; i < npoints1; i++) {
    
    R_CheckUserInterrupt();
    
    dmin = hu;
    d2min = hu2;
    jwhich = -1;
    x1i = x1[i];
    y1i = y1[i];
    z1i = z1[i];
    id1i = id1[i];

    /* search backward from previous nearest neighbour */
    if(lastjwhich > 0) {
      for(j = lastjwhich - 1; j >= 0; --j) {
	dz = z2[j] - z1i;
	dz2 = dz * dz;
	if(dz2 > d2min)
	  break;
#ifdef EXCLUDE
	/* do not compare identical points */
	if(id2[j] != id1i) {
#endif
	  dx = x2[j] - x1i;
	  dy = y2[j] - y1i;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = j;
	  }
#ifdef EXCLUDE
	}
#endif
      }
    }

    /* search forward from previous nearest neighbour  */
    if(lastjwhich < npoints2) {
      for(j = lastjwhich; j < npoints2; ++j) {
	dz = z2[j] - z1i;
	dz2 = dz * dz;
	if(dz2 > d2min)
	  break;
#ifdef EXCLUDE
	/* do not compare identical points */
	if(id2[j] != id1i) {
#endif
	  dx = x2[j] - x1i;
	  dy = y2[j] - y1i;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = j;
	  }
#ifdef EXCLUDE
	}
#endif
      }
    }
#ifdef DIST
    nnd[i] = sqrt(d2min);
#endif
#ifdef WHICH
    /* convert to R indexing */
    nnwhich[i] = jwhich;
#endif
    lastjwhich = jwhich;
  }
}
