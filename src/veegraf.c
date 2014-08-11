/*
  veegraf.c

  $Revision: 1.2 $  $Date: 2013/05/21 08:11:27 $ 

  Given the edges of a graph, determine all "Vees"
  i.e. triples (i, j, k) where i ~ j and i ~ k. 

*/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"

#undef DEBUGVEE

SEXP graphVees(SEXP nv,  /* number of vertices */
	       SEXP iedge,  /* vectors of indices of ends of each edge */   
	       SEXP jedge)  /* all arguments are integer */
/* Edges should NOT be repeated symmetrically. Indices need not be sorted.  */
{
  int Nv, Ne;
  int *ie, *je;         /* edges */
  int *it, *jt, *kt;    /* vectors of indices of triples */ 
  int Nt, Ntmax;        /* number of triples */

  int Nj;
  int *jj; /* scratch storage */

  int i, j, k, m, mj, mk, Nmore, maxchunk;
  
  /* output */
  SEXP iTout, jTout, kTout, out;
  int *ito, *jto, *kto;
  
  /* =================== Protect R objects from garbage collector ======= */
  PROTECT(nv = AS_INTEGER(nv));
  PROTECT(iedge = AS_INTEGER(iedge));
  PROTECT(jedge = AS_INTEGER(jedge));
  /* That's 3 protected objects */

  /* numbers of vertices and edges */
  Nv = *(INTEGER_POINTER(nv)); 
  Ne = LENGTH(iedge);

  /* input arrays */
  ie = INTEGER_POINTER(iedge);
  je = INTEGER_POINTER(jedge);

  /* initialise storage (with a guess at max size) */
  Ntmax = 3 * Ne;
  it = (int *) R_alloc(Ntmax, sizeof(int));
  jt = (int *) R_alloc(Ntmax, sizeof(int));
  kt = (int *) R_alloc(Ntmax, sizeof(int));
  Nt = 0;

  /* initialise scratch storage */
  jj = (int *) R_alloc(Ne, sizeof(int));

  XOUTERCHUNKLOOP(i, 1, Nv, maxchunk, 8196) {
    R_CheckUserInterrupt();
    XINNERCHUNKLOOP(i, 1, Nv, maxchunk, 8196) {

#ifdef DEBUGVEE
      Rprintf("i=%d ---------- \n", i);
#endif

      /* Find Vee triples with apex 'i' */

      /* First, find all vertices j connected to i */
      Nj = 0;
      for(m = 0; m < Ne; m++) {
	if(ie[m] == i) {
	  jj[Nj] = je[m];
	  Nj++;
	} else if(je[m] == i) {
	  jj[Nj] = ie[m];
	  Nj++;
	}
      }

      /* 
	 save triples (i,j,k) 
      */

#ifdef DEBUGVEE
      Rprintf("Nj = %d\n", Nj);
#endif

      if(Nj > 1) {
#ifdef DEBUGVEE
	Rprintf("i=%d\njj=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif

	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    /* add (i, j, k) to list of triangles */
	    if(Nt >= Ntmax) {
	      /* overflow - allocate more space */
	      Nmore = 2 * Ntmax;
#ifdef DEBUGVEE
	      Rprintf("Doubling space from %d to %d\n", Ntmax, Nmore);
#endif
	      it = (int *) S_realloc((char *) it,
				     Nmore,  Ntmax,
				     sizeof(int));
	      jt = (int *) S_realloc((char *) jt,
				     Nmore,  Ntmax,
				     sizeof(int));
	      kt = (int *) S_realloc((char *) kt,
				     Nmore,  Ntmax,
				     sizeof(int));
	      Ntmax = Nmore;
	    }
	    it[Nt] = i;
	    jt[Nt] = j;
	    kt[Nt] = k;
	    Nt++;
	  } 
	}
      }
    }
  }

  /* allocate space for output */
  PROTECT(iTout = NEW_INTEGER(Nt));
  PROTECT(jTout = NEW_INTEGER(Nt));
  PROTECT(kTout = NEW_INTEGER(Nt));
  PROTECT(out   = NEW_LIST(3));
  /* that's 3+4=7 protected objects */
  
  ito = INTEGER_POINTER(iTout);
  jto = INTEGER_POINTER(jTout);
  kto = INTEGER_POINTER(kTout);
  
  /* copy triplet indices to output vectors */
  for(m = 0; m < Nt; m++) {
    ito[m] = it[m];
    jto[m] = jt[m];
    kto[m] = kt[m];
  }
  
  /* insert output vectors in output list */
  SET_VECTOR_ELT(out, 0, iTout);
  SET_VECTOR_ELT(out, 1, jTout);
  SET_VECTOR_ELT(out, 2, kTout);

  UNPROTECT(7);
  return(out);
}
