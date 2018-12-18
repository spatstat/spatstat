
/* 
   mhloop.h

   This file contains the iteration loop 
   for the Metropolis-Hastings algorithm methas.c 

   It is #included several times in methas.c
   with different #defines for the following variables

   MH_MARKED    whether the simulation is marked
               (= the variable 'marked' is TRUE)

   MH_SINGLE    whether there is a single interaction 
              (as opposed to a hybrid of several interactions)

   MH_TEMPER    whether tempering is applied

   MH_TRACKING  whether to save transition history

   MH_DEBUG     whether to print debug information
   
   MH_SNOOP     whether to run visual debugger

   $Revision: 1.23 $  $Date: 2018/12/18 02:43:11 $ 

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#ifndef MH_DEBUG
#define MH_DEBUG NO
#endif

/* ..... Pre-processing: recursively delete illegal/improbable points ..... */

nfree = state.npts - algo.ncond;  /* number of 'free' points */

if(thinstart && nfree > 0) {
  nsuspect = nfree;
  while(nsuspect > 0) {
    /* scan for illegal points */
    ix = state.npts - nsuspect;
    deathprop.ix = ix;
    deathprop.u  = state.x[ix];
    deathprop.v  = state.y[ix];
#if MH_MARKED
    deathprop.mrk = state.marks[ix];
#endif
#if MH_DEBUG
#if MH_MARKED
    Rprintf("check legality of point %d = (%lf, %lf) with mark %d\n", 
	    ix, deathprop.u, deathprop.v, deathprop.mrk);
#else
    Rprintf("check legality of point %d = (%lf, %lf)\n", 
	    ix, deathprop.u, deathprop.v);
#endif
#endif
    /* evaluate conditional intensity without trend terms */

#if MH_SINGLE
    adenom = (*(thecif.eval))(deathprop, state, thecdata);
#else
    adenom = 1.0;
    for(k = 0; k < Ncif; k++)
      adenom *= (*(cif[k].eval))(deathprop, state, cdata[k]);
#endif
#if MH_TEMPER
    adenom = pow(adenom, invtemp);
#endif
#if MH_DEBUG
    Rprintf("cif = %lf\n", adenom);
#endif
    /* accept/reject */
    if(unif_rand() >= adenom) {
#if MH_DEBUG
      Rprintf("deleting illegal/improbable point\n");
#endif
      /* delete point x[ix], y[ix] */
      if(mustupdate) {
	/* Update auxiliary variables first */
#if MH_SINGLE
	(*(thecif.update))(state, deathprop, thecdata);
#else
	for(k = 0; k < Ncif; k++) {
	  if(needupd[k])
	    (*(cif[k].update))(state, deathprop, cdata[k]);
	}
#endif
      }
      state.npts--;
      nfree--;
#if MH_DEBUG
      Rprintf("deleting point %d\n", ix);
      Rprintf("\tnpts=%d\n", state.npts);
#endif
      if(ix < state.npts) {
	for(j = ix; j < state.npts; j++) {
	  state.x[j] = state.x[j+1];
	  state.y[j] = state.y[j+1];
#if MH_MARKED
	  state.marks[j] = state.marks[j+1];
#endif
	}
      }
    }
    nsuspect--;
  }
 }


/* ............... MAIN ITERATION LOOP  ............................. */


OUTERCHUNKLOOP(irep, algo.nrep, maxchunk, 1024) {
  R_CheckUserInterrupt();
  INNERCHUNKLOOP(irep, algo.nrep, maxchunk, 1024) {

#if MH_DEBUG
    Rprintf("iteration %d\n", irep);
#endif

    if(verb) {
      /* print progress message every nverb iterations */
      iverb = irep + 1 + algo.nrep0;
      if((iverb % algo.nverb) == 0)
	Rprintf("iteration %d\n", iverb);
    }

    itype = REJECT;

    nfree = state.npts - algo.ncond;  /* number of 'free' points */

    /* ................  generate proposal ..................... */
    /* Shift or birth/death: */
    if(unif_rand() > algo.p) {
#if MH_DEBUG
      Rprintf("propose birth or death\n");
#endif
      /* Birth/death: */
      if(unif_rand() > algo.q) {
	/* Propose birth: */
	birthprop.u = xpropose[irep];
	birthprop.v = ypropose[irep];
#if MH_MARKED
	birthprop.mrk = mpropose[irep];
#endif
#if MH_DEBUG
#if MH_MARKED
	Rprintf("propose birth at (%lf, %lf) with mark %d\n", 
		birthprop.u, birthprop.v, birthprop.mrk);
#else
	Rprintf("propose birth at (%lf, %lf)\n", birthprop.u, birthprop.v);
#endif
#endif
	/* evaluate conditional intensity */

#if MH_MARKED
	betavalue = model.beta[birthprop.mrk];
#endif

#if MH_SINGLE
	anumer = betavalue * (*(thecif.eval))(birthprop, state, thecdata);
#else
	anumer = betavalue;
	for(k = 0; k < Ncif; k++)
	  anumer *= (*(cif[k].eval))(birthprop, state, cdata[k]);
#endif
#if MH_TEMPER
        anumer = pow(anumer, invtemp);
#endif

	adenom = qnodds*(nfree+1);

#if MH_DEBUG
	Rprintf("cif = %lf, Hastings ratio = %lf\n", anumer, anumer/adenom);
#endif

	/* accept/reject */
	if(unif_rand() * adenom < anumer) {
#if MH_DEBUG
	  Rprintf("accepted birth\n");
#endif
	  itype = BIRTH;  /* Birth proposal accepted. */
	}
#if MH_SNOOP
	/* visual debug */
	mhsnoop(&snooper, irep, &algo, &state, &birthprop, 
		anumer, adenom, &itype);
#endif
#if MH_TRACKING
	/* save transition history */
	if(irep < history.nmax) {
	  history.n++;
	  history.proptype[irep] = BIRTH;
	  history.accepted[irep] = (itype == REJECT) ? 0 : 1;
#ifdef HISTORY_INCLUDES_RATIO
	  history.numerator[irep] = anumer;
	  history.denominator[irep] = adenom;
#endif
	}
#endif
      } else if(nfree > 0) {
	/* Propose death: */
	ix = floor(nfree * unif_rand());
	if(ix < 0) ix = 0;
	ix = algo.ncond + ix;
	if(ix >= state.npts) ix = state.npts - 1;
	deathprop.ix = ix;
	deathprop.u  = state.x[ix];
	deathprop.v  = state.y[ix];
#if MH_MARKED
	deathprop.mrk = state.marks[ix];
#endif
#if MH_DEBUG
#if MH_MARKED
	Rprintf("propose death of point %d = (%lf, %lf) with mark %d\n", 
		ix, deathprop.u, deathprop.v, deathprop.mrk);
#else
	Rprintf("propose death of point %d = (%lf, %lf)\n", 
		ix, deathprop.u, deathprop.v);
#endif
#endif
	/* evaluate conditional intensity */

#if MH_MARKED
	betavalue = model.beta[deathprop.mrk];
#endif

#if MH_SINGLE
	adenom = betavalue * (*(thecif.eval))(deathprop, state, thecdata);
#else
	adenom = betavalue;
	for(k = 0; k < Ncif; k++)
	  adenom *= (*(cif[k].eval))(deathprop, state, cdata[k]);
#endif
#if MH_TEMPER
        adenom = pow(adenom, invtemp);
#endif

	anumer = qnodds * nfree;
#if MH_DEBUG
	Rprintf("cif = %lf, Hastings ratio = %lf\n", adenom, anumer/adenom);
#endif
	/* accept/reject */
	if(unif_rand() * adenom < anumer) {
#if MH_DEBUG
	  Rprintf("accepted death\n");
#endif
	  itype = DEATH; /* Death proposal accepted. */
	}
#if MH_SNOOP
	/* visual debug */
	mhsnoop(&snooper, irep, &algo, &state, &deathprop, 
		anumer, adenom, &itype);
#endif
#if MH_TRACKING
	/* save transition history */
	if(irep < history.nmax) {
	  history.n++;
	  history.proptype[irep] = DEATH;
	  history.accepted[irep] = (itype == REJECT) ? 0 : 1;
#ifdef HISTORY_INCLUDES_RATIO
	  history.numerator[irep] = anumer;
	  history.denominator[irep] = adenom;
#endif
	}
#endif
      }
    } else if(nfree > 0) {
      /* Propose shift: */
      /* point to be shifted */
      ix = floor(nfree * unif_rand());
      if(ix < 0) ix = 0;
      ix = algo.ncond + ix;
      if(ix >= state.npts) ix = state.npts - 1;
      deathprop.ix = ix;
      deathprop.u  = state.x[ix];
      deathprop.v  = state.y[ix];
#if MH_MARKED
      deathprop.mrk = state.marks[ix];
#endif
      /* where to shift */
      permitted = YES;
      shiftprop.ix = ix;
      shiftprop.u = xpropose[irep]; 
      shiftprop.v = ypropose[irep];
#if MH_MARKED
      shiftprop.mrk = mpropose[irep]; 
      if(algo.fixall) permitted = (shiftprop.mrk == deathprop.mrk);
#endif

#if MH_DEBUG
#if MH_MARKED
      Rprintf("propose shift of point %d = (%lf, %lf)[mark %d] to (%lf, %lf)[mark %d]\n", 
	      ix, deathprop.u, deathprop.v, deathprop.mrk, 
	      shiftprop.u, shiftprop.v, shiftprop.mrk);
#else
      Rprintf("propose shift of point %d = (%lf, %lf) to (%lf, %lf)\n", 
	      ix, deathprop.u, deathprop.v, shiftprop.u, shiftprop.v);
#endif
#endif

      /* evaluate cif in two stages */
      cvn = cvd = 1.0;
      if(permitted) {
#if MH_SINGLE
	cvn = (*(thecif.eval))(shiftprop, state, thecdata);
	if(cvn > 0.0) {
	  cvd = (*(thecif.eval))(deathprop, state, thecdata);
	} else {
	  permitted = NO;
	}
#else
	for(k = 0; k < Ncif; k++) {
	  cvn *= (*(cif[k].eval))(shiftprop, state, cdata[k]);
	  if(cvn > 0.0) {
	    cvd *= (*(cif[k].eval))(deathprop, state, cdata[k]);
	  } else {
	    permitted = NO;
	    break; 
	  }
	}
#endif
      } 

      if(permitted) {
#if MH_MARKED
	cvn *= model.beta[shiftprop.mrk];
	cvd *= model.beta[deathprop.mrk];
#endif
#if MH_TEMPER
	cvn = pow(cvn, invtemp);
	cvd = pow(cvd, invtemp);
#endif

#if MH_DEBUG
	Rprintf("cif[old] = %lf, cif[new] = %lf, Hastings ratio = %lf\n", 
		cvd, cvn, cvn/cvd);
#endif
	/* accept/reject */
	if(unif_rand() * cvd < cvn) {
#if MH_DEBUG
	  Rprintf("accepted shift\n");
#endif
	  itype = SHIFT;          /* Shift proposal accepted . */
	}
      } else {
	cvn = 0.0;
	cvd = 1.0;
#if MH_DEBUG
	Rprintf("Forbidden shift");
#endif
      }

#if MH_SNOOP
	/* visual debug */
	mhsnoop(&snooper, irep, &algo, &state, &shiftprop, 
		cvn, cvd, &itype);
#endif
#if MH_TRACKING
      /* save transition history */
      if(irep < history.nmax) {
	history.n++;
	history.proptype[irep] = SHIFT;
	history.accepted[irep] = (itype == REJECT) ? 0 : 1;
#ifdef HISTORY_INCLUDES_RATIO
	  history.numerator[irep] = cvn;
	  history.denominator[irep] = cvd;
#endif
      }
#endif
    }

    if(itype != REJECT) {
      /* ....... implement the transition ............  */
      if(itype == BIRTH) {      
	/* Birth transition */
	/* add point at (u,v) */
#if MH_DEBUG
#if MH_MARKED
	Rprintf("implementing birth at (%lf, %lf) with mark %d\n", 
		birthprop.u, birthprop.v, birthprop.mrk);
#else
	Rprintf("implementing birth at (%lf, %lf)\n", 
		birthprop.u, birthprop.v);
#endif
#endif
	if(state.npts + 1 > state.npmax) {
#if MH_DEBUG
	  Rprintf("!!!!!!!!!!! storage overflow !!!!!!!!!!!!!!!!!\n");
#endif
	  /* storage overflow; allocate more storage */
	  Nmore = 2 * state.npmax;
	  state.x = (double *) S_realloc((char *) state.x, 
					 Nmore,  state.npmax, 
					 sizeof(double));
	  state.y = (double *) S_realloc((char *) state.y, 
					 Nmore,  state.npmax, 
					 sizeof(double));
#if MH_MARKED
	  state.marks = (int *) S_realloc((char *) state.marks, 
					  Nmore,  state.npmax, 
					  sizeof(int));
#endif
	  state.npmax = Nmore;

	  /* call the initialiser again, to allocate additional space */
#if MH_SINGLE
	  thecdata = (*(thecif.init))(state, model, algo);
#else
	  model.ipar = iparvector;
	  for(k = 0; k < Ncif; k++) {
	    if(k > 0)
	      model.ipar += plength[k-1];
	    cdata[k] = (*(cif[k].init))(state, model, algo);
	  }	
#endif
#if MH_DEBUG
	  Rprintf("........... storage extended .................\n");
#endif
	}
	
	if(mustupdate) {
	  /* Update auxiliary variables first */
#if MH_SINGLE
	  (*(thecif.update))(state, birthprop, thecdata);
#else
	  for(k = 0; k < Ncif; k++) {
	    if(needupd[k])
	      (*(cif[k].update))(state, birthprop, cdata[k]);
	  }
#endif
	}
	/* Now add point */
	state.x[state.npts] = birthprop.u;
	state.y[state.npts] = birthprop.v;
#if MH_MARKED
	state.marks[state.npts] = birthprop.mrk;
#endif
	state.npts     = state.npts + 1;
#if MH_DEBUG
	Rprintf("\tnpts=%d\n", state.npts);
#endif
      } else if(itype==DEATH) { 
	/* Death transition */
	/* delete point x[ix], y[ix] */
	if(mustupdate) {
	  /* Update auxiliary variables first */
#if MH_SINGLE
	  (*(thecif.update))(state, deathprop, thecdata);
#else
	  for(k = 0; k < Ncif; k++) {
	    if(needupd[k])
	      (*(cif[k].update))(state, deathprop, cdata[k]);
	  }
#endif
	}
	ix = deathprop.ix;
	state.npts = state.npts - 1;
#if MH_DEBUG
	Rprintf("implementing death of point %d\n", ix);
	Rprintf("\tnpts=%d\n", state.npts);
#endif
	if(ix < state.npts) {
	  for(j = ix; j < state.npts; j++) {
	    state.x[j] = state.x[j+1];
	    state.y[j] = state.y[j+1];
#if MH_MARKED
	    state.marks[j] = state.marks[j+1];
#endif
	  }
	}
      } else {
	/* Shift transition */
	/* Shift (x[ix], y[ix]) to (u,v) */
#if MH_DEBUG
#if MH_MARKED
	Rprintf("implementing shift from %d = (%lf, %lf)[%d] to (%lf, %lf)[%d]\n", 
		deathprop.ix, deathprop.u, deathprop.v, deathprop.mrk,
		shiftprop.u, shiftprop.v, shiftprop.mrk);
#else
	Rprintf("implementing shift from %d = (%lf, %lf) to (%lf, %lf)\n", 
		deathprop.ix, deathprop.u, deathprop.v,
		shiftprop.u, shiftprop.v);
	Rprintf("\tnpts=%d\n", state.npts);
#endif
#endif
	if(mustupdate) {
	  /* Update auxiliary variables first */
#if MH_SINGLE
	  (*(thecif.update))(state, shiftprop, thecdata);
#else
	  for(k = 0; k < Ncif; k++) {
	    if(needupd[k])
	      (*(cif[k].update))(state, shiftprop, cdata[k]);
	  }
#endif
	}
	ix = shiftprop.ix;
	state.x[ix] = shiftprop.u;
	state.y[ix] = shiftprop.v;
#if MH_MARKED
	state.marks[ix] = shiftprop.mrk;
#endif
      }
#if MH_DEBUG
    } else {
      Rprintf("rejected\n");
#endif
    }
  }
}
