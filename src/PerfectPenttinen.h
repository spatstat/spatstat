
// ........................... Penttinen process ................
// $Revision: 1.2 $  $Date: 2016/02/02 01:30:01 $

class PenttProcess : public PointProcess {
 public:
  double beta, gamma, radius, reachsquared, loggamma2pi;
  int ishard;
  PenttProcess(double xmin, double xmax, double ymin, double ymax, 
	       double b, double g, double r);
  ~PenttProcess(){}
  void NewEvent(double *x, double *y, char *InWindow);
  void GeneratePoisson(Point *headPoint, 
			       long int *GeneratedPoints,
			       long int *LivingPoints,
			       long int *NoP);
  double Interaction(double dsquared);
};

PenttProcess::PenttProcess(double xmin, double xmax, 
			   double ymin, double ymax, 
			   double b, double g, double r) :
  PointProcess(xmin, xmax, ymin, ymax){
    beta = b; gamma = g; radius = r; 
    ishard = (gamma <= DOUBLE_EPS);
    loggamma2pi = M_2PI * (ishard? 0.0 : log(gamma));
    reachsquared = 4.0 * radius * radius;
    InteractionRange = 2.0 * radius;
    TotalBirthRate = beta*(xmax-xmin)*(ymax-ymin);
  }  

double PenttProcess::Interaction(double dsquared)
{
  double rtn, z, z2;
  rtn = 1.0;
  if(dsquared < reachsquared) {
    if(ishard) return(0.0);
    z2 = dsquared/reachsquared;
    z = sqrt(z2);
    if(z < 1.0) {
      rtn = exp(loggamma2pi * (acos(z) - z * sqrt(1.0 - z2)));
    }
  }
  return(rtn);
}

void PenttProcess::NewEvent(double *x, double *y, char *InWindow)
{
  double Xdim, Ydim;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  *x = slumptal()*Xdim+Xmin;
  *y = slumptal()*Ydim+Ymin;
  *InWindow = 1;
}

void PenttProcess::GeneratePoisson(Point *headPoint, 
			      long int *GeneratedPoints,
			      long int *LivingPoints,
			      long int *NoP)
{
  int i;
  double xtemp, ytemp, L, Xdim, Ydim;
  struct Point *TempPoint;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  L = beta*Xdim*Ydim;
  *GeneratedPoints = poisson(L);
  *LivingPoints = *GeneratedPoints;
  for (i=1; i<=*GeneratedPoints ; i++){
    //Rprintf("Generating PenttProcess Poisson 3\n");
    //scanf("%f",&f1);
    xtemp = slumptal()*Xdim+Xmin;
    ytemp = slumptal()*Ydim+Ymin;
    //
    //Rprintf("Generating PenttProcess Poisson 3.2\n");
    TempPoint = ALLOCATE(struct Point);
    //
    TempPoint->X = xtemp;
    TempPoint->Y = ytemp;
    TempPoint->No = i;
    TempPoint->R = slumptal();
    //Rprintf("Generating PenttProcess Poisson 3.6\n");
    TempPoint->next = headPoint->next;
    headPoint->next = TempPoint;
    *NoP = *NoP + 1;
  }
}

// ........................... Interface to R ..........................

extern "C" {
  SEXP PerfectPenttinen(SEXP beta,
			SEXP gamma,
			SEXP r,
			SEXP xrange,
			SEXP yrange) {

    // input parameters
    double Beta, Gamma, R, Xmin, Xmax, Ymin, Ymax;
    double *Xrange, *Yrange;
    // internal
    int xcells, ycells;
    long int StartTime, EndTime;
    // output 
    int noutmax;
    SEXP xout, yout, nout, out;
    double *xx, *yy;
    int *nn;

    // protect arguments from garbage collector    
    PROTECT(beta   = AS_NUMERIC(beta));
    PROTECT(gamma   = AS_NUMERIC(gamma));
    PROTECT(r      = AS_NUMERIC(r));
    PROTECT(xrange = AS_NUMERIC(xrange));
    PROTECT(yrange = AS_NUMERIC(yrange));
    // that's 5 protected objects

    // extract arguments
    Beta   = *(NUMERIC_POINTER(beta));
    Gamma  = *(NUMERIC_POINTER(gamma));
    R      = *(NUMERIC_POINTER(r));

    Xrange = NUMERIC_POINTER(xrange);
    Xmin   = Xrange[0];
    Xmax   = Xrange[1];
    Yrange = NUMERIC_POINTER(yrange);
    Ymin   = Yrange[0];
    Ymax   = Yrange[1];

    // compute cell array size
    xcells = (int) floor((Xmax-Xmin)/ R);
    if(xcells > 9) xcells = 9; if(xcells < 1) xcells = 1;
    ycells = (int) floor((Ymax-Ymin)/ R);

    Xrange = NUMERIC_POINTER(xrange);
    Xmin   = Xrange[0];
    Xmax   = Xrange[1];
    Yrange = NUMERIC_POINTER(yrange);
    Ymin   = Yrange[0];
    Ymax   = Yrange[1];

    // compute cell array size
    xcells = (int) floor((Xmax-Xmin)/ R);
    if(xcells > 9) xcells = 9; if(xcells < 1) xcells = 1;
    ycells = (int) floor((Ymax-Ymin)/ R);
    if(ycells > 9) ycells = 9; if(ycells < 1) ycells = 1;
#ifdef DBGS
    Rprintf("xcells %d   ycells %d\n",xcells,ycells);
    Rprintf("Initialising\n");
#endif

    // Initialise Penttinen point process
    PenttProcess ExampleProcess(Xmin,Xmax,Ymin,Ymax,Beta,Gamma,R);  
    // Initialise point pattern
    Point2Pattern ExamplePattern(Xmin,Xmax,Ymin,Ymax, xcells, ycells);
    // parameters: min x, max x, min y, max y, "cells" in x and y direction
    // used for speeding up neighbour counting, 9 is max here
    
#ifdef DBGS
    Rprintf("Initialisation complete\n");
#endif

    // Synchronise random number generator 
    GetRNGstate();

    // Initialise perfect sampler
    Sampler PerfectSampler(&ExampleProcess);
    
    // Perform perfect sampling
    PerfectSampler.Sim(&ExamplePattern, &StartTime, &EndTime);
    
    // Synchronise random number generator 
    PutRNGstate();

    // Get upper estimate of number of points
    noutmax = ExamplePattern.UpperCount() + 1;
    
    // Allocate space for output
    PROTECT(xout = NEW_NUMERIC(noutmax));
    PROTECT(yout = NEW_NUMERIC(noutmax));
    PROTECT(nout = NEW_INTEGER(1));
    xx = NUMERIC_POINTER(xout);
    yy = NUMERIC_POINTER(yout);
    nn = INTEGER_POINTER(nout);

    // copy data into output storage
    ExamplePattern.Return(xx, yy, nn, noutmax);

    // pack up into output list
    PROTECT(out  = NEW_LIST(3));
    SET_VECTOR_ELT(out, 0, xout);
    SET_VECTOR_ELT(out, 1, yout);
    SET_VECTOR_ELT(out, 2, nout);
    
    // return 
    UNPROTECT(9);  // 5 arguments plus xout, yout, nout, out
    return(out);
  }
}
