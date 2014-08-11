#include <R.h>

/* 
   linnncross.c

   Shortest-path distances between nearest neighbours in linear network
   One pattern to another pattern

   $Revision: 1.1 $  $Date: 2013/10/21 02:01:29 $

   linndcross      
   linndxcross     

*/

#define DPATH(I,J) dpath[(I) + Nv * (J)]
#define ANSWER(I,J) answer[(I) + Np * (J)]
#define EUCLID(X,Y,U,V) sqrt(pow((X)-(U),2)+pow((Y)-(V),2))

/* definition of linndcross */
#define FNAME linndcross
#undef  EXCLU
#define WHICH

#include "linnncross.h"

#undef  FNAME
#undef  EXCLU
#undef  WHICH

/* definition of linndxcross */

#define FNAME linndxcross
#define EXCLU
#define WHICH

#include "linnncross.h"
