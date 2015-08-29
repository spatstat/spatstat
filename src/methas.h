/* 
   Definitions of types and data structures for Metropolis-Hastings

   State       Current state of point pattern

   Model       Model parameters passed from R

   Cdata       (pointer to) model parameters and precomputed data in C

   Algor       Algorithm parameters (p, q, nrep etc)

   Propo       Proposal in Metropolis-Hastings algorithm

   History     Transition history of MH algorithm

   Cifns       Set of functions for computing the conditional intensity
               for a point process model. 
	       This consists of three functions
                    init(State, Model, Algor) .... initialises auxiliary data
		    eval(State, Propo) ........... evaluates cif
		    update(State,Propo) .......... updates auxiliary data

 */


/* Current state of point pattern */
typedef struct State { 
  double *x;     /* vectors of Cartesian coordinates */
  double *y;
  int *marks;    /* vector of mark values */
  int npts;       /* current number of points */
  int npmax;      /* storage limit */
  int ismarked;   /* whether the pattern is marked */
} State;

/* Parameters of model passed from R */
typedef struct Model {
  double *beta;     /* vector of activity parameters */
  double *ipar;     /* vector of interaction parameters */
  double *period;  /* width & height of rectangle, if torus */
  int ntypes;      /* number of possible marks */
} Model;

/* 
   A pointer to Cdata 
   is a pointer to C storage for parameters of model
*/

typedef void Cdata;

/* RMH Algorithm parameters */
typedef struct Algor {
  double p;         /* probability of proposing shift */
  double q;         /* conditional probability of proposing death */
  int fixall;       /* if TRUE, only shifts of location are feasible */
  int ncond;        /* For conditional simulation, 
		       the first 'ncond' points are fixed */
  int nrep;        /* number of iterations */
  int nverb;       /* print report every 'nverb' iterations */
  int nrep0;       /* number of iterations already performed 
		      in previous blocks - for reporting purposes */
  int tempered;    /* TRUE if tempering is applied */
  double invtemp;  /* inverse temperature if tempering is applied */
} Algor;

/* Metropolis-Hastings proposal */
typedef struct Propo {
  double u;         /* location of point of interest */
  double v;
  int mrk;       /* mark of point of interest */
  int ix;           /* index of point of interest, if already in pattern */
  int itype;        /* transition type */
} Propo;

/* transition codes 'itype' */
#define REJECT 0
#define BIRTH 1
#define DEATH 2
#define SHIFT 3

#define HISTORY_INCLUDES_RATIO

/* Record of transition history */
typedef struct History {
  int nmax;              /* length of vectors */
  int n;                 /* number of events recorded */
  int *proptype;         /* vector: proposal type */
  int *accepted;         /* vector: 0 for reject, 1 for accept */
#ifdef HISTORY_INCLUDES_RATIO
  double *numerator;     /* vectors: Hastings ratio numerator & denominator  */
  double *denominator;
#endif
} History;

/* conditional intensity functions */

typedef Cdata *   (*initfunptr)(State state, Model model, Algor algo);
typedef double    (*evalfunptr)(Propo prop,  State state, Cdata *cdata);
typedef void      (*updafunptr)(State state, Propo prop,  Cdata *cdata);

typedef struct Cifns {
  initfunptr init;
  evalfunptr eval;
  updafunptr update;
  int        marked;
} Cifns;

#define NEED_UPDATE(X) ((X).update != (updafunptr) NULL)

#define NULL_CIFNS { (initfunptr) NULL, (evalfunptr) NULL, (updafunptr) NULL, NO}

/* miscellaneous macros */

#include "yesno.h"

# define MAT(X,I,J,M) (X[(I)+(J)*(M)])




