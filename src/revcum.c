/* 

   revcum.c

   $Revision: 1.3 $  $Date: 2016/12/30 07:28:13 $

   Reverse cumulative sums

*/

void drevcumsum(double *x, int *nx) {
  int i;
  double sumx;
  double *xp;
  
  i = *nx - 1;
  xp = x + i;
  sumx = *xp;
  while(i > 0) {
    --i;
    --xp;
    sumx += *xp;
    *xp = sumx;
  }
}

void irevcumsum(int *x, int *nx) {
  int i;
  int sumx;
  int *xp;
  
  i = *nx - 1;
  xp = x + i;
  sumx = *xp;
  while(i > 0) {
    --i;
    --xp;
    sumx += *xp;
    *xp = sumx;
  }
}
