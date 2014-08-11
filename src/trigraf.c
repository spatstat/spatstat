/*

  trigraf.c

  Form list of all triangles in a planar graph, given list of edges
  
  $Revision: 1.13 $     $Date: 2012/04/06 09:26:50 $

  Form list of all triangles in a planar graph, given list of edges

  Note: vertex indices ie, je are indices in R.
        They are handled without converting to C convention,
        because we only need to test equality and ordering.
	(*except in 'trioxgraph'*)

  Called by .C:
  -------------
  trigraf()  Generic C implementation with fixed storage limit
             usable with Delaunay triangulation

  trigrafS() Faster version when input data are sorted
	     (again with predetermined storage limit)
	     suited for handling Delaunay triangulation

  Called by .Call:
  ---------------
  trigraph()   Version with dynamic storage allocation

  triograph()  Faster version assuming 'iedge' is sorted in increasing order

  trioxgraph()  Even faster version for use with quadrature schemes

  Diameters:
  -----------
  triDgraph() Also computes diameters of triangles

*/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"

#undef DEBUGTRI

void trigraf(nv, ne, ie, je, ntmax, nt, it, jt, kt, status)
     /* inputs */
     int *nv;         /* number of graph vertices */
     int *ne;         /* number of edges */
     int *ie, *je;    /* vectors of indices of ends of each edge */ 
     int *ntmax;      /* length of storage space for triangles */
     /* output */
     int *nt;              /* number of triangles (<= *ntmax) */
     int *it, *jt, *kt;    /* vectors of indices of vertices of triangles */ 
     int *status;          /* 0 if OK, 1 if overflow */
{
  int Nv, Ne, Ntmax;
  int Nt, Nj, m, i, j, k, mj, mk, maxchunk;
  int *jj;
  
  Nv = *nv;
  Ne = *ne;
  Ntmax = *ntmax;

  /* initialise scratch storage */
  jj = (int *) R_alloc(Ne, sizeof(int));
  Nt = 0;

  /* vertex index i ranges from 1 to Nv */
  XOUTERCHUNKLOOP(i, 1, Nv, maxchunk, 8196) {
    R_CheckUserInterrupt();
    XINNERCHUNKLOOP(i, 1, Nv, maxchunk, 8196) {
      /* Find triangles involving vertex 'i'
	 in which 'i' is the lowest-numbered vertex */

      /* First, find vertices j > i connected to i */
      Nj = 0;
      for(m = 0; m < Ne; m++) {
	if(ie[m] == i) {
	  j = je[m];
	  if(j > i) {
	    jj[Nj] = j;
	    Nj++;
	  }
	} else if(je[m] == i) {
	  j = ie[m];
	  if(j > i) {
	    jj[Nj] = j;
	    Nj++;
	  }
	}
      }

      /* 
	 Determine which pairs of vertices j, k are joined by an edge;
	 save triangles (i,j,k) 
      */

      if(Nj > 1) {
	/* Sort jj in ascending order */
	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(k < j) {
	      /* swap */
	      jj[mk] = j;
	      jj[mj] = k;
	      j = k;
	    }
	  }
	}
	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(j != k) {
	      /* Run through edges to determine whether j, k are neighbours */
	      for(m = 0; m < Ne; m++) {
		if((ie[m] == j && je[m] == k)
		   || (ie[m] == k && je[m] == j)) {
		  /* add (i, j, k) to list of triangles */
		  if(Nt >= Ntmax) {
		    /* overflow - exit */
		    *status = 1;
		    return;
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
      }
    }
  }
  *nt = Nt;
  *status = 0;
}


/* faster version of trigraf() 
   assuming that 
            ie[m] < je[m]
            ie[] is in ascending order
            je[] is in ascending order within ie[],
	          that is, je[ie[]=i] is in ascending order for each fixed i
*/

void trigrafS(nv, ne, ie, je, ntmax, nt, it, jt, kt, status)
     /* inputs */
     int *nv;         /* number of graph vertices */
     int *ne;         /* number of edges */
     int *ie, *je;    /* vectors of indices of ends of each edge */ 
     int *ntmax;      /* length of storage space for triangles */
     /* output */
     int *nt;              /* number of triangles */
     int *it, *jt, *kt;    /* vectors of indices of vertices of triangles */ 
     int *status;          /* 0 if OK, 1 if overflow */
{
  int Ne, Nt, Ntmax;
  int m, i, j, k, mj, mk;
  int firstedge, lastedge;
  
  Ne = *ne;
  Ntmax = *ntmax;

  /* nv is not used, but retained for harmony with trigraf */
  /* Avoid compiler warnings */
  Nt = *nv;

  /* initialise output */
  Nt = 0;

  lastedge = -1;
  while(lastedge + 1 < Ne) {
    if(lastedge % 256 == 0) R_CheckUserInterrupt();
    /* 
       Consider next vertex i.
       The edges (i,j) with i < j appear contiguously in the edge list.
    */
    firstedge = lastedge + 1;
    i = ie[firstedge]; 
    for(m= firstedge+1; m < Ne && ie[m] == i; m++)
      ;
    lastedge = m-1;
    /* 
       Consider each pair j, k of neighbours of i, where i < j < k. 
       Scan entire edge list to determine whether j, k are joined by an edge.
       If so, save triangle (i,j,k) 
    */
    if(lastedge > firstedge) {
      for(mj = firstedge; mj < lastedge; mj++) {
	j = je[mj];
	for(mk = firstedge+1; mk <= lastedge; mk++) {
	  k = je[mk];
	  /* Run through edges to determine whether j, k are neighbours */
	  for(m = 0; m < Ne && ie[m] < j; m++) 
	    ;
	  while(m < Ne && ie[m] == j) {
	    if(je[m] == k) {
	      /* add (i, j, k) to list of triangles */
	      if(Nt >= Ntmax) {
		/* overflow - exit */
		*status = 1;
		return;
	      }
	      it[Nt] = i;
	      jt[Nt] = j;
	      kt[Nt] = k;
	      Nt++;
	    }
	    m++;
	  }
	}
      }
    }
  }

  *nt = Nt;
  *status = 0;
}


/* ------------------- callable by .Call ------------------------- */


SEXP trigraph(SEXP nv,  /* number of vertices */
	      SEXP iedge,  /* vectors of indices of ends of each edge */   
	      SEXP jedge)  /* all arguments are integer */
{
  int Nv, Ne;
  int *ie, *je;         /* edges */
  int *it, *jt, *kt;    /* vectors of indices of vertices of triangles */ 
  int Nt, Ntmax;        /* number of triangles */

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

#ifdef DEBUGTRI
      Rprintf("i=%d ---------- \n", i);
#endif

      /* Find triangles involving vertex 'i'
	 in which 'i' is the lowest-numbered vertex */

      /* First, find vertices j > i connected to i */
      Nj = 0;
      for(m = 0; m < Ne; m++) {
	if(ie[m] == i) {
	  j = je[m];
	  if(j > i) {
	    jj[Nj] = j;
	    Nj++;
	  }
	} else if(je[m] == i) {
	  j = ie[m];
	  if(j > i) {
	    jj[Nj] = j;
	    Nj++;
	  }
	}
      }

      /* 
	 Determine which pairs of vertices j, k are joined by an edge;
	 save triangles (i,j,k) 
      */

#ifdef DEBUGTRI
      Rprintf("Nj = %d\n", Nj);
#endif

      if(Nj > 1) {
#ifdef DEBUGTRI
	Rprintf("i=%d\njj=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif
	/* Sort jj in ascending order */
	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(k < j) {
	      /* swap */
	      jj[mk] = j;
	      jj[mj] = k;
	      j = k;
	    }
	  }
	}
#ifdef DEBUGTRI
	Rprintf("sorted=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif

	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(j != k) {
	      /* Run through edges to determine whether j, k are neighbours */
	      for(m = 0; m < Ne; m++) {
		if((ie[m] == j && je[m] == k)
		   || (ie[m] == k && je[m] == j)) {
		  /* add (i, j, k) to list of triangles */
		  if(Nt >= Ntmax) {
		    /* overflow - allocate more space */
		    Nmore = 2 * Ntmax;
#ifdef DEBUGTRI
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
		  /* output indices in R convention */
		  it[Nt] = i;
		  jt[Nt] = j;
		  kt[Nt] = k;
		  Nt++;
		}
	      }
	    }
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
  
  /* copy triangle indices to output vectors */
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


/* faster version assuming iedge is in increasing order */

SEXP triograph(SEXP nv,  /* number of vertices */
	       SEXP iedge,  /* vectors of indices of ends of each edge */   
	       SEXP jedge)  /* all arguments are integer */
{
  int Nv, Ne;
  int *ie, *je;         /* edges */
  int *it, *jt, *kt;    /* vectors of indices of vertices of triangles */ 
  int Nt, Ntmax;        /* number of triangles */

  int Nj;
  int *jj; /* scratch storage */

  int i, j, k, m, mj, mk, maxjk, Nmore, maxchunk;
  
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

#ifdef DEBUGTRI
      Rprintf("i=%d ---------- \n", i);
#endif

      /* Find triangles involving vertex 'i'
	 in which 'i' is the lowest-numbered vertex */

      /* First, find vertices j > i connected to i */
      Nj = 0;
      for(m = 0; m < Ne; m++) {
	if(ie[m] == i) {
	  j = je[m];
	  if(j > i) {
	    jj[Nj] = j;
	    Nj++;
	  }
	} else if(je[m] == i) {
	  j = ie[m];
	  if(j > i) {
	    jj[Nj] = j;
	    Nj++;
	  }
	}
      }

      /* 
	 Determine which pairs of vertices j, k are joined by an edge;
	 save triangles (i,j,k) 
      */

#ifdef DEBUGTRI
      Rprintf("Nj = %d\n", Nj);
#endif

      if(Nj > 1) {
#ifdef DEBUGTRI
	Rprintf("i=%d\njj=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif
	/* Sort jj in ascending order */
	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(k < j) {
	      /* swap */
	      jj[mk] = j;
	      jj[mj] = k;
	      j = k;
	    }
	  }
	}
#ifdef DEBUGTRI
	Rprintf("sorted=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif

	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(j != k) {
	      /* Run through edges to determine whether j, k are neighbours */
	      maxjk = (j > k) ? j : k;
	      for(m = 0; m < Ne; m++) {
		if(ie[m] > maxjk) break;
		/* 
		   since iedge is in increasing order, the test below
		   will always be FALSE when ie[m] > max(j,k)
		*/
		if((ie[m] == j && je[m] == k)
		   || (ie[m] == k && je[m] == j)) {
		  /* add (i, j, k) to list of triangles */
		  if(Nt >= Ntmax) {
		    /* overflow - allocate more space */
		    Nmore = 2 * Ntmax;
#ifdef DEBUGTRI
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
  
  /* copy triangle indices to output vectors */
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

/* 
   Even faster version using information about dummy vertices.
   Dummy-to-dummy edges are forbidden.

   For generic purposes use 'friendly' for 'isdata'
   Edge between j and k is possible iff friendly[j] || friendly[k].
   Edges with friendly = FALSE cannot be connected to one another.

 */


SEXP trioxgraph(SEXP nv,  /* number of vertices */
		SEXP iedge,  /* vectors of indices of ends of each edge */   
		SEXP jedge,
		SEXP friendly)  /* indicator vector, length nv */
{
  /* input */
  int Nv, Ne;
  int *ie, *je;         /* edges */
  int *friend;         /* indicator */

  /* output */
  int *it, *jt, *kt;    /* vectors of indices of vertices of triangles */ 
  int Nt, Ntmax;        /* number of triangles */

  /* scratch storage */
  int Nj;
  int *jj; 
  int i, j, k, m, mj, mk, maxjk, Nmore, maxchunk;
  
  /* output to R */
  SEXP iTout, jTout, kTout, out;
  int *ito, *jto, *kto;
  
  /* =================== Protect R objects from garbage collector ======= */
  PROTECT(nv = AS_INTEGER(nv));
  PROTECT(iedge = AS_INTEGER(iedge));
  PROTECT(jedge = AS_INTEGER(jedge));
  PROTECT(friendly = AS_INTEGER(friendly));
  /* That's 4 protected objects */

  /* numbers of vertices and edges */
  Nv = *(INTEGER_POINTER(nv)); 
  Ne = LENGTH(iedge);

  /* input arrays */
  ie = INTEGER_POINTER(iedge);
  je = INTEGER_POINTER(jedge);
  friend = INTEGER_POINTER(friendly);

  /* initialise storage (with a guess at max size) */
  Ntmax = 3 * Ne;
  it = (int *) R_alloc(Ntmax, sizeof(int));
  jt = (int *) R_alloc(Ntmax, sizeof(int));
  kt = (int *) R_alloc(Ntmax, sizeof(int));
  Nt = 0;

  /* initialise scratch storage */
  jj = (int *) R_alloc(Ne, sizeof(int));

  /* convert to C indexing convention */
  for(m = 0; m < Ne; m++) {
    ie[m] -= 1;
    je[m] -= 1;
  }

  OUTERCHUNKLOOP(i, Nv, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Nv, maxchunk, 8196) {

#ifdef DEBUGTRI
      Rprintf("i=%d ---------- \n", i);
#endif

      /* Find triangles involving vertex 'i'
	 in which 'i' is the lowest-numbered vertex */

      /* First, find vertices j > i connected to i */
      Nj = 0;
      for(m = 0; m < Ne; m++) {
	if(ie[m] == i) {
	  j = je[m];
	  if(j > i) {
	    jj[Nj] = j;
	    Nj++;
	  }
	} else if(je[m] == i) {
	  j = ie[m];
	  if(j > i) {
	    jj[Nj] = j;
	    Nj++;
	  }
	}
      }

      /* 
	 Determine which pairs of vertices j, k are joined by an edge;
	 save triangles (i,j,k) 
      */

#ifdef DEBUGTRI
      Rprintf("Nj = %d\n", Nj);
#endif

      if(Nj > 1) {
#ifdef DEBUGTRI
	Rprintf("i=%d\njj=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif
	/* Sort jj in ascending order */
	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(k < j) {
	      /* swap */
	      jj[mk] = j;
	      jj[mj] = k;
	      j = k;
	    }
	  }
	}
#ifdef DEBUGTRI
	Rprintf("sorted=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif

	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(j != k && (friend[j] || friend[k])) {
	      /* Run through edges to determine whether j, k are neighbours */
	      maxjk = (j > k) ? j : k;
	      for(m = 0; m < Ne; m++) {
		if(ie[m] > maxjk) break;
		/* 
		   since iedge is in increasing order, the test below
		   will always be FALSE when ie[m] > max(j,k)
		*/
		if((ie[m] == j && je[m] == k)
		   || (ie[m] == k && je[m] == j)) {
		  /* add (i, j, k) to list of triangles */
		  if(Nt >= Ntmax) {
		    /* overflow - allocate more space */
		    Nmore = 2 * Ntmax;
#ifdef DEBUGTRI
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
		  /* convert back to R indexing */
		  it[Nt] = i + 1;
		  jt[Nt] = j + 1;
		  kt[Nt] = k + 1;
		  Nt++;
		} 
	      }
	    }
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
  /* that's 4+4=8 protected objects */
  
  ito = INTEGER_POINTER(iTout);
  jto = INTEGER_POINTER(jTout);
  kto = INTEGER_POINTER(kTout);
  
  /* copy triangle indices to output vectors */
  for(m = 0; m < Nt; m++) {
    ito[m] = it[m];
    jto[m] = jt[m];
    kto[m] = kt[m];
  }
  
  /* insert output vectors in output list */
  SET_VECTOR_ELT(out, 0, iTout);
  SET_VECTOR_ELT(out, 1, jTout);
  SET_VECTOR_ELT(out, 2, kTout);

  UNPROTECT(8);
  return(out);
}

/* 
   also calculates diameter (max edge length) of triangle
*/

SEXP triDgraph(SEXP nv,  /* number of vertices */
	       SEXP iedge,  /* vectors of indices of ends of each edge */   
	       SEXP jedge,
	       SEXP edgelength)   /* edge lengths */
{
  int Nv, Ne;
  int *ie, *je;         /* edges */
  double *edgelen;      

  int *it, *jt, *kt;    /* vectors of indices of vertices of triangles */ 
  double *dt;           /* diameters (max edge lengths) of triangles */
  int Nt, Ntmax;        /* number of triangles */

  /* scratch storage */
  int Nj;
  int *jj; 
  double *dd;

  int i, j, k, m, mj, mk, Nmore, maxchunk;
  double dij, dik, djk, diam;
  
  /* output */
  SEXP iTout, jTout, kTout, dTout, out;
  int *ito, *jto, *kto;
  double *dto;
  
  /* =================== Protect R objects from garbage collector ======= */
  PROTECT(nv = AS_INTEGER(nv));
  PROTECT(iedge = AS_INTEGER(iedge));
  PROTECT(jedge = AS_INTEGER(jedge));
  PROTECT(edgelength = AS_NUMERIC(edgelength));
  /* That's 4 protected objects */

  /* numbers of vertices and edges */
  Nv = *(INTEGER_POINTER(nv)); 
  Ne = LENGTH(iedge);

  /* input arrays */
  ie = INTEGER_POINTER(iedge);
  je = INTEGER_POINTER(jedge);
  edgelen = NUMERIC_POINTER(edgelength);

  /* initialise storage (with a guess at max size) */
  Ntmax = 3 * Ne;
  it = (int *) R_alloc(Ntmax, sizeof(int));
  jt = (int *) R_alloc(Ntmax, sizeof(int));
  kt = (int *) R_alloc(Ntmax, sizeof(int));
  dt = (double *) R_alloc(Ntmax, sizeof(double));
  Nt = 0;

  /* initialise scratch storage */
  jj = (int *) R_alloc(Ne, sizeof(int));
  dd = (double *) R_alloc(Ne, sizeof(double));

  XOUTERCHUNKLOOP(i, 1, Nv, maxchunk, 8196) {
    R_CheckUserInterrupt();
    XINNERCHUNKLOOP(i, 1, Nv, maxchunk, 8196) {

#ifdef DEBUGTRI
      Rprintf("i=%d ---------- \n", i);
#endif

      /* Find triangles involving vertex 'i'
	 in which 'i' is the lowest-numbered vertex */

      /* First, find vertices j > i connected to i */
      Nj = 0;
      for(m = 0; m < Ne; m++) {
	if(ie[m] == i) {
	  j = je[m];
	  if(j > i) {
	    jj[Nj] = j;
	    dd[Nj] = edgelen[m];
	    Nj++;
	  }
	} else if(je[m] == i) {
	  j = ie[m];
	  if(j > i) {
	    jj[Nj] = j;
	    dd[Nj] = edgelen[m];
	    Nj++;
	  }
	}
      }

      /* 
	 Determine which pairs of vertices j, k are joined by an edge;
	 save triangles (i,j,k) 
      */

#ifdef DEBUGTRI
      Rprintf("Nj = %d\n", Nj);
#endif

      if(Nj > 1) {
#ifdef DEBUGTRI
	Rprintf("i=%d\njj=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif
	/* Sort jj in ascending order */
	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    if(k < j) {
	      /* swap */
	      jj[mk] = j;
	      jj[mj] = k;
	      dik = dd[mj];
	      dd[mj] = dd[mk];
	      dd[mk] = dik;
	      j = k;
	    }
	  }
	}
#ifdef DEBUGTRI
	Rprintf("sorted=\n", i);
	for(mj = 0; mj < Nj; mj++) Rprintf("%d ", jj[mj]);
	Rprintf("\n\n");
#endif

	for(mj = 0; mj < Nj-1; mj++) {
	  j = jj[mj];
	  dij = dd[mj];
	  for(mk = mj+1; mk < Nj; mk++) {
	    k = jj[mk];
	    dik = dd[mk];
	    if(j != k) {
	      /* Run through edges to determine whether j, k are neighbours */
	      for(m = 0; m < Ne; m++) {
		if((ie[m] == j && je[m] == k)
		   || (ie[m] == k && je[m] == j)) {
		  /* triangle (i, j, k) */
		  /* determine triangle diameter */
		  diam = (dij > dik) ? dij : dik;
		  djk = edgelen[m];
		  if(djk > diam) diam = djk; 
		  /* add (i, j, k) to list of triangles */
		  if(Nt >= Ntmax) {
		    /* overflow - allocate more space */
		    Nmore = 2 * Ntmax;
#ifdef DEBUGTRI
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
		    dt = (double *) S_realloc((char *) dt,
					      Nmore,  Ntmax,
					      sizeof(double));
		    Ntmax = Nmore;
		  }
		  it[Nt] = i;
		  jt[Nt] = j;
		  kt[Nt] = k;
		  dt[Nt] = diam; 
		  Nt++;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /* allocate space for output */
  PROTECT(iTout = NEW_INTEGER(Nt));
  PROTECT(jTout = NEW_INTEGER(Nt));
  PROTECT(kTout = NEW_INTEGER(Nt));
  PROTECT(dTout = NEW_NUMERIC(Nt));
  PROTECT(out   = NEW_LIST(4));
  /* that's 4+5=9 protected objects */
  
  ito = INTEGER_POINTER(iTout);
  jto = INTEGER_POINTER(jTout);
  kto = INTEGER_POINTER(kTout);
  dto = NUMERIC_POINTER(dTout);
  
  /* copy triangle indices to output vectors */
  for(m = 0; m < Nt; m++) {
    ito[m] = it[m];
    jto[m] = jt[m];
    kto[m] = kt[m];
    dto[m] = dt[m];
  }
  
  /* insert output vectors in output list */
  SET_VECTOR_ELT(out, 0, iTout);
  SET_VECTOR_ELT(out, 1, jTout);
  SET_VECTOR_ELT(out, 2, kTout);
  SET_VECTOR_ELT(out, 3, dTout);

  UNPROTECT(9);
  return(out);
}

