#include <R.h>
#include <Rdefines.h>

/*

  fastinterv.c

  Fast version of findInterval when breaks are known to be evenly spaced
  and are known to embrace the data.

  $Revision: 1.2 $ $Date: 2015/10/19 11:09:18 $

*/

void fastinterv(x, n, brange, nintervals, y) 
     double *x, *brange;
     int *n, *nintervals;
     int *y;
{
  double bmin, bmax, db;
  int i, j, m, N;

  m = *nintervals;
  N = *n;

  bmin = brange[0];
  bmax = brange[1];
  db = (bmax - bmin)/m;
  
  for(i = 0; i < N; i++) {
    j = (int) ceil((x[i] - bmin)/db);
    if(j <= 0) { j = 1; } else if(j > m) { j = m; }
    y[i] = j;
  }
}

	       
