/*
  quasirandom.c 

  Quasi-random sequence generators

  Copyright (C) Adrian Baddeley 2014
  GNU Public Licence version 2 | 3

  $Revision: 1.1 $  $Date: 2014/03/17 03:31:59 $

*/

#include <math.h>

void Corput(base, n, result) 
     int *base, *n;
     double *result; 
{
  int b, N, i, j;
  register double f, f0, z;

  N = *n;
  b = *base;

  f0 = 1.0/((double) b);

  for(i = 0; i < N; i++) {
    j = i+1;
    z = 0;
    f = f0;
    while(j > 0) {
      z = z + f * (j % b);
      j = j/b; 
      f = f / ((double) b);
    }
    result[i] = z;
  }
}
