/*

   dinfty.c

   $Revision: 1.6 $   $Date: 2011/09/20 07:42:18 $

   Code by Dominic Schuhmacher

   Modified by Adrian Baddeley

*/

#include <stdio.h>
#include <R.h>

#define COST(I,J) (d)[n * (J) + (I)]

int arraymax(int *a, int n);
void swap(int i, int j, int *a);
int largestmobpos(int *mobile, int *current, int *collectvals, int n);

/* ------------ The main function ----------------------------- */

void dinfty_R(int *d, int *num, int *assignment) {
   int i,j; /* indices */
   int lmp, lmq; /* largest mobile position and its neighbor */
   int newmax;
   int n, currmin;
   int *current, *travel, *mobile, *assig, *distrelev, *collectvals;

   n = *num;

   /* scratch space */
   assig = (int *) R_alloc((long) n, sizeof(int)); 
   travel = (int *) R_alloc((long) n, sizeof(int)); 
   mobile = (int *) R_alloc((long) n, sizeof(int)); 
   current = (int *) R_alloc((long) n, sizeof(int)); 
   distrelev = (int *) R_alloc((long) n, sizeof(int));

   collectvals = (int *) R_alloc((long) (n * n), sizeof(int));


/*                                                               */
/* We use the Johnson-Trotter Algorithm for listing permutations */
/*                                                               */

/* Initialize the algorithm */
   for (i = 0; i < n; i++) {
      travel[i] = -1;   /* all numbers traveling to the left */
      mobile[i] = 1;    /* all numbers mobile */
      current[i] = i;   /* current permutation is the identity */
      assig[i] = i;     /* best permutation up to now is the identity */
      distrelev[i] = COST(i, i);   /* pick relevant entries in the cost matrix */
   }
   currmin = arraymax(distrelev, n);   /* minimal max up to now */

/* The main loop */
   while(arraymax(mobile, n) == 1) {
     lmp = largestmobpos(mobile, current, collectvals, n);
      lmq = lmp + travel[lmp];
      swap(lmp, lmq, current);
      swap(lmp, lmq, travel);
      for (i = 0; i < n; i++) {
         if (current[i] > current[lmq])
            travel[i] = -travel[i];
         j = i + travel[i];
         if (j < 0 || j > n-1 || current[i] < current[j])
            mobile[i] = 0;
         else
            mobile[i] = 1;
         distrelev[i] = COST(i, current[i]);
      }
      /* Calculation of new maximal value */
      newmax = arraymax(distrelev, n);
      if (newmax < currmin) {
         currmin = newmax;
         for (i = 0; i < n; i++) {
            assig[i] = current[i];
         }
      }
   }
/* For testing: print distance from within C program
   Rprintf("Prohorov distance is %d\n", currmin);     */

/* "Return" the final assignment */
   for (i = 0; i < n; i++) {
      assignment[i] = assig[i] + 1;
   }

}


/* ------------------------------------------------------------*/


/* Maximal element of an integer array */
int arraymax(int *a, int n) {
  int i, amax;
  if(n < 1)
    return(-1);
  amax = a[0];
  if(n > 1)
    for(i = 0; i < n; i++)
      if(a[i] > amax) amax = a[i];
  return(amax);
}


/* Swap elements i and j in array a */

void swap(int i, int j, int *a) {
   int v;

   v = a[i];
   a[i] = a[j];
   a[j] = v;
}


/* Return index of largest mobile number in current */
int largestmobpos(int *mobile, int *current, int *collectvals, int n) {
   int i,j, maxval;

   j = 0;
   for (i = 0; i < n; i++) {
      if (mobile[i] == 1) {
         collectvals[j] = current[i];
         j++;
      }
   }
   maxval = arraymax(collectvals, j);
   for (i = 0; i < n; i++) {
      if (current[i] == maxval) {
         return(i);
      }
   }
   error("Internal error: largestmobpos failed");
   return(0);
}
