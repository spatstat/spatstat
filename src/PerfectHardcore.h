
// ........................... Hardcore process ..........................
// $Revision: 1.4 $  $Date: 2012/03/10 11:23:09 $

class HardcoreProcess : public PointProcess {
 public:
  double beta, R, Rsquared;
  HardcoreProcess(double xmin, double xmax, double ymin, double ymax, 
		double b, double Ri);
  ~HardcoreProcess(){}
  void NewEvent(double *x, double *y, char *InWindow);
  void GeneratePoisson(Point *headPoint, 
			       long int *GeneratedPoints,
			       long int *LivingPoints,
			       long int *NoP);
  double Interaction(double dsquared);
};

HardcoreProcess::HardcoreProcess(double xmin, double xmax, 
			      double ymin, double ymax, 
			      double b, double Ri) :
  PointProcess(xmin, xmax, ymin, ymax){
  beta = b; R = Ri; 
  Rsquared = R * R; 
  InteractionRange = R;
  TotalBirthRate = beta*(xmax-xmin)*(ymax-ymin);
  }  

double HardcoreProcess::Interaction(double dsquared)
{
  double rtn;
  rtn = 1;
  if(dsquared < Rsquared) rtn = 0;
  return(rtn);
}

void HardcoreProcess::NewEvent(double *x, double *y, char *InWindow)
{
  double Xdim, Ydim;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  *x = slumptal()*Xdim+Xmin;
  *y = slumptal()*Ydim+Ymin;
  *InWindow = 1;
}

void HardcoreProcess::GeneratePoisson(Point *headPoint, 
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
    //Rprintf("Generating HardcoreProcess Poisson 3\n");
    //scanf("%f",&f1);
    xtemp = slumptal()*Xdim+Xmin;
    ytemp = slumptal()*Ydim+Ymin;
    //
    //Rprintf("Generating HardcoreProcess Poisson 3.2\n");
    TempPoint = ALLOCATE(struct Point);
    //
    TempPoint->X = xtemp;
    TempPoint->Y = ytemp;
    TempPoint->No = i;
    TempPoint->R = slumptal();
    //Rprintf("Generating HardcoreProcess Poisson 3.6\n");
    TempPoint->next = headPoint->next;
    headPoint->next = TempPoint;
    *NoP = *NoP + 1;
  }
}

// ........................... Interface to R ..........................

extern "C" {
  SEXP PerfectHardcore(SEXP beta,
		      SEXP r,
		      SEXP xrange,
		      SEXP yrange) {

    // input parameters
    double Beta, R, Xmin, Xmax, Ymin, Ymax;
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
    PROTECT(r      = AS_NUMERIC(r));
    PROTECT(xrange = AS_NUMERIC(xrange));
    PROTECT(yrange = AS_NUMERIC(yrange));
    // that's 4 protected objects

    // extract arguments
    Beta   = *(NUMERIC_POINTER(beta));
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
    if(ycells > 9) ycells = 9; if(ycells < 1) ycells = 1;
#ifdef DBGS
    Rprintf("xcells %d   ycells %d\n",xcells,ycells);
    Rprintf("Initialising\n");
#endif

    // Initialise Hardcore point process
    HardcoreProcess ExampleProcess(Xmin,Xmax,Ymin,Ymax, Beta, R);  
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
    UNPROTECT(8);  // 4 arguments plus xout, yout, nout, out
    return(out);
  }
}
