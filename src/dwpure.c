/*

   dwpure.c

   $Revision: 1.5 $   $Date: 2011/09/20 07:54:53 $

   Code by Dominic Schuhmacher
   
*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>

typedef struct State {
  int n1, n2; 
  /* vectors of length n1 (rows) and n2 (cols) */
  int *rowmass, *colmass;  /* mass to be moved from row / to col */
  int *rowlab, *collab;  /* row and col labels
         (specify previous node (row for collab, col for rowlab)) */
  int *rowflow, *colflow;  /* second component of labels 
         (specify flow through current node) */
  int *rowsurplus, *colsurplus;  /* the surplus in each row/col under the current flow */
  int *dualu, *dualv;     /* vectors of dual variables (u for rows, v for cols) */
  int *rowhelper, *colhelper;   /* helping vector to store intermediate results */
      /* could be local in initcost at the moment */
  /* n by n matrices */
  int *d;                /* matrix of costs */ 
  int *flowmatrix;        /* matrix of flows */
  int *arcmatrix;  /* matrix of arcs for restriced primal problem 
         (1 if arc, 0 if no arc) should be unsigned char to save memory
          however need to workout problem with R_alloc first (see below) */
   /* n*n vector */
   int *collectvals;
} State;

#define COST(I,J,STATE,NVALUE) ((STATE)->d)[(NVALUE) * (J) + (I)]
#define FLOW(I,J,STATE,NVALUE) ((STATE)->flowmatrix)[(NVALUE) * (J) + (I)]
#define ARC(I,J,STATE,NVALUE) ((STATE)->arcmatrix)[(NVALUE) * (J) + (I)]
#define MIN(A,B) ((A)<(B) ? (A) : (B))

int arraysum(int *a, int n);
int arraymin(int *a, int n);
void initvalues(State *state);
void maxflow(State *state);
void updateduals(State *state);
void augmentflow(int startcol, State *state);

/* ------------ The main function ----------------------------- */

void dwpure(int *d, int *rmass, int *cmass, int *numr, int *numc, int *flowmatrix)
{
   int i,j; /* indices */
   int n1,n2;
   unsigned char feasible = 0; /* boolean for main loop */
   State state;

   /* inputs */
   state.n1 = n1 = *numr;
   state.n2 = n2 = *numc;
   state.d = d;
   state.rowmass = rmass;
   state.colmass = cmass;
   /* scratch space */
   state.rowlab = (int *) R_alloc((long) n1, sizeof(int));
   state.collab = (int *) R_alloc((long) n2, sizeof(int));
   state.rowflow = (int *) R_alloc((long) n1, sizeof(int));
   state.colflow = (int *) R_alloc((long) n2, sizeof(int));
   state.rowsurplus = (int *) R_alloc((long) n1, sizeof(int));
   state.colsurplus = (int *) R_alloc((long) n2, sizeof(int));
   state.dualu = (int *) R_alloc((long) n1, sizeof(int));
   state.dualv = (int *) R_alloc((long) n2, sizeof(int));
   state.rowhelper = (int *) R_alloc((long) n1, sizeof(int));
   state.colhelper = (int *) R_alloc((long) n2, sizeof(int));
   state.flowmatrix = (int *) R_alloc((long) (n1 * n2), sizeof(int));
   state.arcmatrix = (int *) R_alloc((long) (n1 * n2), sizeof(int));
   state.collectvals = (int *) R_alloc((long) (n1 * n2), sizeof(int));

   for (i = 0; i < n1; ++i) {
   for (j = 0; j < n2; ++j) {
      state.flowmatrix[(n1)*(j) + i] = 0;
      state.arcmatrix[(n1)*(j) + i] = 0;
      state.collectvals[(n1)*(j) + i] = 0;
   }
   }
   for (i = 0; i < n1; ++i) {
      state.rowlab[i] = 0;
      state.rowflow[i] = 0;
      state.rowsurplus[i] = 0;
      state.dualu[i] = 0;
      state.rowhelper[i] = 0;
   }
   for (j = 0; j < n2; ++j) {
      state.collab[j] = 0;
      state.colflow[j] = 0;
      state.colsurplus[j] = 0;
      state.dualv[j] = 0;
      state.colhelper[j] = 0;
   }


/* Initialize dual variables, arcmatrix, and surpluses */
   initvalues(&state);

/* For testing: print out cost matrix 
   for (i = 0; i < n1; ++i) {
   for (j = 0; j < n2; ++j) {
      Rprintf("%d ", COST(i, j, &state, n1));
   }
   Rprintf("\n");
   }   */

/* The main loop */
   while(feasible == 0) {
      maxflow(&state);
      if (arraysum(state.rowsurplus, n1) > 0) {
         updateduals(&state);  /* also updates arcmatrix */
      }
      else {
         feasible = 1;
      }
   }

/* "Return" the final flowmatrix */
   for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
      flowmatrix[n1*j+i] = state.flowmatrix[n1*j+i];
      }
   }
}


/* ------------ Functions called by dwpure_R ------------------------- */


/* Sum of integer array */
int arraysum(int *a, int n) {
   int i;
   int asum = 0;
   for (i = 0; i < n; i++)
      asum += a[i];
   return(asum);
}

/* Minimal element of an integer array */
int arraymin(int *a, int n) {
  int i, amin;
  if (n < 1)
    return(-1);
  amin = a[0];
  if (n > 1)
    for (i = 0; i < n; i++)
      if (a[i] < amin) amin = a[i];
  return(amin);
}


/* Initialize cost matrix: subtract in each row its minimal entry (from all the
entries in the row), then subtract in each column its minimal entry (from all the
entries in the column) */
void initvalues(State *state) {
   int i,j,n1,n2;

   n1 = state->n1;
   n2 = state->n2;

   /* Initial surpluses; can I do this shorter? later on surpluses are updated in
      flow augmentation step */
   for (i = 0; i < n1; i++)
      state->rowsurplus[i] = state->rowmass[i];
   for (j = 0; j < n2; j++)
      state->colsurplus[j] = state->colmass[j];

   for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) 
         state->colhelper[j] = COST(i, j, state, n1);
      state->dualu[i] = arraymin(state->colhelper, n2);
   }
   for (j = 0; j < n2; j++) {
      for (i = 0; i < n1; i++) 
	 state->rowhelper[i] = COST(i, j, state, n1) - state->dualu[i];
      state->dualv[j] = arraymin(state->rowhelper, n1);
   }
   for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
         if (COST(i, j, state, n1) == state->dualu[i] + state->dualv[j])
            ARC(i, j, state, n1) = 1;
         else
            ARC(i, j, state, n1) = 0;
      }
   }
}

/* Maximize the flow on the (zeros of the) current cost matrix */
void maxflow(State *state) {
   int breakthrough; /* col. no. in which breakthrough occurs */
   unsigned char labelfound = 1; /* 0 if no more labels can be found */
   int i,j,n1,n2;

   n1 = state->n1;
   n2 = state->n2;

   while (labelfound == 1) {
      breakthrough = -1;
      /* initialize labels */
      for (i = 0; i < n1; i++) {
         if (state->rowsurplus[i] > 0) {
            state->rowlab[i] = -5;
            state->rowflow[i] = state->rowsurplus[i];
         }
         else {
            state->rowlab[i] = -1;  /* setting rowflow to zero isn't necessary! */
         }
      }
      for (j = 0; j < n2; j++)
         state->collab[j] = -1;   /* setting colflow to zero isn't necessary! */
      /* -1 means "no index", -5 means "source label" (rows only) */

      while (labelfound == 1 && breakthrough == -1) {
         labelfound = 0;
         /* label unlabeled column j that permits flow from some labeled row i */
         /* ("permits flow" means arcmatrix[i][j] = 1). Do so for every j */
         for (i = 0; i < n1; i++) {
            if (state->rowlab[i] != -1) {
               for (j = 0; j < n2; j++) {
                  if (ARC(i, j, state, n1) == 1 && state->collab[j] == -1) {
                     state->collab[j] = i;
                     state->colflow[j] = state->rowflow[i];
                     labelfound = 1;
                     if (state->colsurplus[j] > 0 && breakthrough == -1)
                        breakthrough = j;
                  }
               }
            }
         }
         /* label unlabeled row i that already sends flow to some labeled col j */
         /* ("already sends" means flowmatrix[i][j] > 0). Do so for every i             */
         for (j = 0; j < n2; j++) {
            if (state->collab[j] != -1) {
               for (i = 0; i < n1; i++) {
                  if (FLOW(i, j, state, n1) > 0 && state->rowlab[i] == -1) {
                     state->rowlab[i] = j;
                     state->rowflow[i] = MIN(state->colflow[j],FLOW(i, j, state, n1));
                     labelfound = 1;
                  }
               }
            }
         }
      }
      if (breakthrough != -1) augmentflow(breakthrough, state);
   }
}


/* Update the dual variables (called if solution of restricted primal is not feasible
for the original problem): determine the minimum over the submatrix given by all
labeled rows and unlabeled columns, and subtract it from all labeled rows and add
it to all labeled columns. */
void updateduals(State *state) 
{
   int i,j,n1,n2,mini;
   int count = 0; 

   n1 = state->n1;
   n2 = state->n2;

   for (i = 0; i < n1; i++) {
     for (j = 0; j < n2; j++) {
       if (state->rowlab[i] != -1 && state->collab[j] == -1) {
	 state->collectvals[count] = COST(i, j, state, n1) - state->dualu[i] - state->dualv[j];
	 count++;
       }
     }
   }
   mini = arraymin(state->collectvals, count);
   for (i = 0; i < n1; i++) {
     if (state->rowlab[i] != -1)
       state->dualu[i] += mini;
   }
   for (j = 0; j < n2; j++){
     if (state->collab[j] != -1)
       state->dualv[j] -= mini;
   }
   for (i = 0; i < n1; i++) {
     for (j = 0; j < n2; j++) {
       if (COST(i, j, state, n1) == state->dualu[i] + state->dualv[j])
         ARC(i, j, state, n1) = 1;
       else
         ARC(i, j, state, n1) = 0;
     }
   }

}

/* Augment the flow on the graph given by arcmatrix (by aug)
according to the row and column labels starting in column startcol */
/* Adjust the surpluses while we're at it (first row and last col have -aug) */ 
void augmentflow(int startcol, State *state) {
   int k,l,aug,n1;
  /* int i,j,k,l,aug,n1,n2; */

   n1 = state->n1;

   l = startcol;
   aug = MIN(state->colflow[l], state->colsurplus[l]);
   state->colsurplus[l] -= aug;

   k = state->collab[l];
   FLOW(k, l, state, n1) += aug;
   l = state->rowlab[k];
   while (l != -5) {
      FLOW(k, l, state, n1) -= aug;
      k = state->collab[l];
      FLOW(k, l, state, n1) += aug;
      l = state->rowlab[k];
   }
   state->rowsurplus[k] -= aug;
}
