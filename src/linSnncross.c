#include <R.h>
#include "yesno.h"

/* 
   linSnncross.c

   Shortest-path distances between nearest neighbours in linear network
   One pattern to another pattern

   $Revision: 1.4 $  $Date: 2018/12/18 02:43:11 $

   'Sparse version' 

   Works with sparse representation
   Does not allow 'exclusion'
   Requires point data to be ordered by segment index.

   linSnndcross      
   linSnndwhich

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

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



