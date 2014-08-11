
// ........................... Strauss process ..........................
//  $Revision: 1.3 $ $Date: 2012/03/10 11:52:48 $

class StraussProcess : public PointProcess {
 public:
  double beta, gamma, R, Rsquared;
  StraussProcess(double xmin, double xmax, double ymin, double ymax, 
		double b, double g, double Ri);
  ~StraussProcess(){}
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

StraussProcess::StraussProcess(double xmin, double xmax, 
			      double ymin, double ymax, 
			      double b, double g, double Ri) :
  PointProcess(xmin, xmax, ymin, ymax){
  beta = b; gamma = g; R = Ri; 
  Rsquared = R * R; 
  InteractionRange = R;
  TotalBirthRate = beta*(xmax-xmin)*(ymax-ymin);
  }  

double StraussProcess::Interaction(double dsquared)
{
  double rtn;
  rtn = 1;
  if(dsquared < Rsquared) rtn = gamma;
  return(rtn);
}

void StraussProcess::NewEvent(double *x, double *y, char *InWindow)
{
  double Xdim, Ydim;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  *x = slumptal()*Xdim+Xmin;
  *y = slumptal()*Ydim+Ymin;
  *InWindow = 1;
}

void StraussProcess::GeneratePoisson(Point *headPoint, 
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
    //Rprintf("Generating StraussProcess Poisson 3\n");
    //scanf("%f",&f1);
    xtemp = slumptal()*Xdim+Xmin;
    ytemp = slumptal()*Ydim+Ymin;
    //
    //Rprintf("Generating StraussProcess Poisson 3.2\n");
    TempPoint = ALLOCATE(struct Point);
    //
    TempPoint->X = xtemp;
    TempPoint->Y = ytemp;
    TempPoint->No = i;
    TempPoint->R = slumptal();
    //Rprintf("Generating StraussProcess Poisson 3.6\n");
    TempPoint->next = headPoint->next;
    headPoint->next = TempPoint;
    *NoP = *NoP + 1;
  }
}

//void StraussProcess::CalcBeta(long int xsidepomm, long int ysidepomm, 
//		   double *betapomm){ 
//  long int i,j,k;
//  k=0;
//  //  Rprintf("\ndiagnostic message: Strauss CalcBeta... %ld %ld\n",xsidepomm,ysidepomm);
//  for(i=0; i<xsidepomm; i++){
//    for(j=0; j<ysidepomm; j++){
//      *(betapomm + i*ysidepomm + j) = this->beta;
//      k++;
//    }
//  } 
//}

//void StraussProcess::CheckBeta(long int xsidepomm, long int ysidepomm, 
//		   double *betapomm){ 
//  long int i,j,k;
//  //  double d1;
//  k=0;
//  //  Rprintf("\ndiagnostic message: Strauss CalcBeta... %ld %ld\n",xsidepomm,ysidepomm);
//  for(i=0; i<xsidepomm; i++){
//    for(j=0; j<ysidepomm; j++){
//      if((fabs(*(betapomm + i*ysidepomm + j)- beta)>0.001) && (k==0)){
//	Rprintf("%f %f %f %ld %ld\n",fabs(*(betapomm + i*ysidepomm + j)- beta),
//	       *(betapomm + i*ysidepomm + j),beta,i,j);
//	k++;
//	//	scanf("%lf",&d1);
//      }
//    }
//  } 
//}

//double StraussProcess::lnCondInt(struct Point2 *TempCell, 
//				 Point2Pattern *p2p){
//  double f1;
//  long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx,k;
//  double dy,dx, lnCI,dst2;
//  struct Point2 *TempCell2;
//  
//  f1 = (TempCell->X-p2p->Xmin)/p2p->XCellDim;  xc = int(f1);
//  CLAMP(xc, 0, p2p->MaxXCell, "xc");
//  f1 = (TempCell->Y-p2p->Ymin)/p2p->YCellDim;  yc = int(f1);
//  CLAMP(yc, 0, p2p->MaxYCell, "yc");
//  
//  dx = (Xmax-Xmin)/(double(p2p->MaxXCell+1));
//  dy = (Ymax-Ymin)/(double(p2p->MaxYCell+1));
//  rx = int(this->InteractionRange/dx+1.0);
//  ry = int(this->InteractionRange/dy+1.0);
//  
//  lnCI = log(TempCell->Beta);
//
//  k = 0;
//  
//  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
//  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
//  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
//  if((yc-ry)>=0) fy=yc-ry; else fy = 0;
//
//  //Rprintf("MCI! %d %d %d %d\n",fx,tx,fy,ty);
//
//  for(xco = fx; xco <= tx; xco++){
//    for(yco = fy; yco <= ty; yco++){
//      CHECK(p2p->headCell[xco][yco], 
//	    "internal error: p2p->headCell[xco][yco] is null in lnCondInt()");
//      TempCell2 = p2p->headCell[xco][yco]->next;
//      CHECK(TempCell2, "internal error: TempCell2 is null in lnCondInt()");
//      while(TempCell2!=TempCell2->next){
//	if(TempCell2 != TempCell){
//	  k++;
//	  dst2 = pow(TempCell->X-TempCell2->X,2)+
//	        pow(TempCell->Y-TempCell2->Y,2);
//	  lnCI += log(Interaction(dst2));
//	}
//	TempCell2 = TempCell2->next; 
//	CHECK(TempCell2, 
//	      "internal error: TempCell2 is null in lnCondInt() loop");
//      }
//    }
//  }
//  return(lnCI);
//}

//void StraussProcess::Beta(struct Point2 *TempCell){
//  TempCell->Beta = beta;
//}

//void StraussProcess::CalcBeta(Point2Pattern *p2p){
//  long int xco,yco;
//  //  double dy,dx;
//  struct Point2 *TempMother;
//
//  for(xco = 0; xco <= p2p->MaxXCell; xco++){
//    for(yco = 0; yco <= p2p->MaxYCell; yco++){
//      CHECK(p2p->headCell[xco][yco], 
//	    "internal error: p2p->headCell[xco][yco] is null in CalcBeta()");
//      TempMother = p2p->headCell[xco][yco]->next;
//      CHECK(TempMother, "internal error: TempMother is null in CalcBeta()");
//      while(TempMother!=TempMother->next){
//	TempMother->Beta = this->beta;
//	TempMother = TempMother->next;
//	CHECK(TempMother, 
//	      "internal error: TempMother is null in CalcBeta() loop");
//      }
//    }
//  }
//}

// ........................... Interface to R ..........................

extern "C" {
  SEXP PerfectStrauss(SEXP beta,
		      SEXP gamma,
		      SEXP r,
		      SEXP xrange,
		      SEXP yrange) {

    // input parameters
    double Beta, Gamma, R, Xmin, Xmax, Ymin, Ymax;
    double *Xrange, *Yrange;
    // internal
    int xcells, ycells;
    PointProcess *TheProcess;
    long int EndTime, StartTime;
    // output 
    int noutmax;
    SEXP xout, yout, nout, out;
    double *xx, *yy;
    int *nn;

    SEXP stout, etout;
    int *ss, *ee;

    // protect arguments from garbage collector    
    PROTECT(beta   = AS_NUMERIC(beta));
    PROTECT(gamma  = AS_NUMERIC(gamma));
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
    if(ycells > 9) ycells = 9; if(ycells < 1) ycells = 1;
#ifdef DBGS
    Rprintf("xcells %d   ycells %d\n",xcells,ycells);
    Rprintf("Initialising\n");
#endif

    // Initialise Strauss point process
    StraussProcess ExampleProcess(Xmin,Xmax,Ymin,Ymax, Beta, Gamma, R);  

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
    PROTECT(stout = NEW_INTEGER(1));
    PROTECT(etout = NEW_INTEGER(1));
    xx = NUMERIC_POINTER(xout);
    yy = NUMERIC_POINTER(yout);
    nn = INTEGER_POINTER(nout);
    ss = INTEGER_POINTER(stout);
    ee = INTEGER_POINTER(etout);

    // copy data into output storage
    ExamplePattern.Return(xx, yy, nn, noutmax);
    *ss = StartTime;
    *ee = EndTime;

    // pack up into output list
    PROTECT(out  = NEW_LIST(5));
    SET_VECTOR_ELT(out, 0, xout);
    SET_VECTOR_ELT(out, 1, yout);
    SET_VECTOR_ELT(out, 2, nout);
    SET_VECTOR_ELT(out, 3, stout);
    SET_VECTOR_ELT(out, 4, etout);
    
    // return 
    UNPROTECT(11);  // 5 arguments plus xout, yout, nout, stout, etout, out
    return(out);
  }
}
