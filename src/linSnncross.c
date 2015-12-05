#include <R.h>
#include "yesno.h"

/* 
   linSnncross.c

   Shortest-path distances between nearest neighbours in linear network
   One pattern to another pattern

   $Revision: 1.3 $  $Date: 2015/11/28 10:08:55 $

   'Sparse version' 

   Works with sparse representation
   Does not allow 'exclusion'
   Requires point data to be ordered by segment index.

   linSnndcross      
   linSnndwhich

*/

void Clinvdist(), Clinvwhichdist();  /* functions from linvdist.c */

#undef HUH

/* definition of linSnndcross */
#define FNAME linSnndcross
#undef WHICH
#include "linSnncross.h"

/* definition of linSnndwhich */
#undef  FNAME
#define FNAME linSnndwhich
#define WHICH
#include "linSnncross.h"



