
// ..................... Strauss-Hardcore process ..........................
//  $Revision: 1.2 $ $Date: 2012/03/10 11:23:17 $

class StraussHardProcess : public PointProcess {
 public:
  double beta, gamma, H, R, Hsquared, Rsquared;
  StraussHardProcess(double xmin, double xmax, double ymin, double ymax, 
		     double b, double g, double Ri, double Hc);
  ~StraussHardProcess(){}
  void NewEvent(double *x, double *y, char *InWindow);
  void GeneratePoisson(Point *headPoint, 
			       long int *GeneratedPoints,
			       long int *LivingPoints,
			       long int *NoP);
  double Interaction(double dsquared);
  //  void CalcBeta(long int xsidepomm, long int ysidepomm, 
  //	   double *betapomm);
  //  void CheckBeta(long int xsidepomm, long int ysidepomm, 
  //		 double *betapomm);
  //  double lnCondInt(struct Point2 *TempCell, Point2Pattern *p2p);
  //  void Beta(struct Point2 *TempCell);
  //  void CalcBeta(Point2Pattern *p2p);
};

StraussHardProcess::StraussHardProcess(double xmin, double xmax, 
			      double ymin, double ymax, 
			      double b, double g, double Ri, double Hc) :
  PointProcess(xmin, xmax, ymin, ymax){
    beta = b; gamma = g; R = Ri;  H = Hc; 
    Rsquared = R * R; 
    Hsquared = H * H; 
    InteractionRange = R;
    TotalBirthRate = beta*(xmax-xmin)*(ymax-ymin);
  }  

double StraussHardProcess::Interaction(double dsquared)
{
  if(dsquared >= Rsquared) return(1.0);
  if(dsquared >= Hsquared) return(gamma);
  return(0.0);
}

void StraussHardProcess::NewEvent(double *x, double *y, char *InWindow)
{
  double Xdim, Ydim;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  *x = slumptal()*Xdim+Xmin;
  *y = slumptal()*Ydim+Ymin;
  *InWindow = 1;
}

void StraussHardProcess::GeneratePoisson(Point *headPoint, 
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
    //Rprintf("Generating StraussHardProcess Poisson 3\n");
    //scanf("%f",&f1);
    xtemp = slumptal()*Xdim+Xmin;
    ytemp = slumptal()*Ydim+Ymin;
    //
    //Rprintf("Generating StraussHardProcess Poisson 3.2\n");
    TempPoint = ALLOCATE(struct Point);
    //
    TempPoint->X = xtemp;
    TempPoint->Y = ytemp;
    TempPoint->No = i;
    TempPoint->R = slumptal();
    //Rprintf("Generating StraussHardProcess Poisson 3.6\n");
    TempPoint->next = headPoint->next;
    headPoint->next = TempPoint;
    *NoP = *NoP + 1;
  }
}

// ........................... Interface to R ..........................

extern "C" {
  SEXP PerfectStraussHard(SEXP beta,
		      SEXP gamma,
		      SEXP r,
		      SEXP hc,
		      SEXP xrange,
		      SEXP yrange) {

    // input parameters
    double Beta, Gamma, R, H, Xmin, Xmax, Ymin, Ymax;
    double *Xrange, *Yrange;
    // internal
    int xcells, ycells;
    PointProcess *TheProcess;
    long int StartTime, EndTime;
    // output 
    int noutmax;
    SEXP xout, yout, nout, out;
    double *xx, *yy;
    int *nn;

    // protect arguments from garbage collector    
    PROTECT(beta   = AS_NUMERIC(beta));
    PROTECT(gamma  = AS_NUMERIC(gamma));
    PROTECT(r      = AS_NUMERIC(r));
    PROTECT(hc     = AS_NUMERIC(hc));
    PROTECT(xrange = AS_NUMERIC(xrange));
    PROTECT(yrange = AS_NUMERIC(yrange));
    // that's 6 protected objects

    // extract arguments
    Beta   = *(NUMERIC_POINTER(beta));
    Gamma  = *(NUMERIC_POINTER(gamma));
    R      = *(NUMERIC_POINTER(r));
    H      = *(NUMERIC_POINTER(hc));
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

    // Initialise StraussHard point process
    StraussHardProcess ExampleProcess(Xmin,Xmax,Ymin,Ymax, Beta, Gamma, R, H);  

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
    UNPROTECT(10);  // 6 arguments plus xout, yout, nout, out
    return(out);
  }
}
