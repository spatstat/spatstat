//  Debug switch 
//  #define DBGS 

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R_ext/Constants.h>

// #include <stdio.h>
// FILE *out;
// File i/o is deprecated in R implementation


#ifdef DBGS
#define CHECK(PTR,MESSAGE) if(((void *) PTR) == ((void *) NULL)) error(MESSAGE)
#define CLAMP(X, LOW, HIGH, XNAME) \
  if((X) > (HIGH)) { \
     Rprintf("Value of %s exceeds upper limit %d\n", XNAME, HIGH); \
     X = HIGH; \
  } else if((X) < (LOW)) { \
     Rprintf("Value of %s is below %d\n", XNAME, LOW); \
     X = LOW; \
  } 
#else
#define CHECK(PTR,MESSAGE)
#define CLAMP(X, LOW, HIGH, XNAME) \
  if((X) > (HIGH)) X = HIGH; else if((X) < (LOW)) X = LOW; 
#endif

// .........................................
// memory allocation 
// using R_alloc

#define ALLOCATE(TYPE)  (TYPE *) R_alloc(1, sizeof(TYPE))
#define FREE(PTR) 

// Alternative using Calloc and Free
// #define ALLOCATE(TYPE)  (TYPE *) Calloc(1, sizeof(TYPE))
// #define FREE(PTR) Free(PTR)

void R_CheckUserInterrupt(void);

struct Point{ long int No; float X; float Y; float R; struct Point *next; }; 

struct Point2{ long int No; float X; float Y; 
  char InLower[2]; 
  double Beta; double TempBeta; struct Point2 *next; }; 

struct Point3{ char Case; char XCell; char YCell; struct Point3 *next; }; 

// const float Pi=3.141593;

double slumptal(void){
  return(runif((double) 0.0, (double) 1.0));
}

long int poisson(double lambda){
  return((long int)rpois(lambda));
}

// ........................... Point patterns ..........................

class Point2Pattern {
public:
  long int UpperLiving[2];
  long int MaxXCell, MaxYCell, NoP;
  double XCellDim, YCellDim, Xmin, Xmax, Ymin, Ymax;
  struct Point2 *headCell[10][10],*dummyCell;
  char DirX[10], DirY[10];
 
  Point2Pattern(double xmin, double xmax,
		double ymin, double ymax, 
		long int mxc, long int myc){
    long int i,j;
    UpperLiving[0] = 0;
    UpperLiving[1] = 0;
    Xmin = xmin; Xmax = xmax;
    Ymin = ymin; Ymax = ymax;
    DirX[1] = 1; DirY[1] = 0;
    DirX[2] = 1; DirY[2] = -1;
    DirX[3] = 0; DirY[3] = -1;
    DirX[4] = -1; DirY[4] = -1;
    DirX[5] = -1; DirY[5] = 0;
    DirX[6] = -1; DirY[6] = 1;
    DirX[7] = 0; DirY[7] = 1;
    DirX[8] = 1; DirY[8] = 1;    
    NoP = 0;
    //
    dummyCell = ALLOCATE(struct Point2);
    //
    dummyCell->next = dummyCell;
    dummyCell->No = 0;
    MaxXCell = mxc; MaxYCell = myc;
    if(MaxXCell>9) MaxXCell = 9;
    if(MaxYCell>9) MaxYCell = 9;
    for(i=0;i<=MaxXCell;i++){
      for(j=0;j<=MaxYCell;j++){
	//
	headCell[i][j] = ALLOCATE(struct Point2);
	//
	headCell[i][j]->next=dummyCell;
      }
    }
    XCellDim = (Xmax-Xmin)/((double)(MaxXCell+1));
    YCellDim = (Ymax-Ymin)/((double)(MaxYCell+1));
  };
  ~Point2Pattern(){}
  //  void Print();
  void Return(double *X, double *Y, int *num, int maxnum);
  long int Count();
  long int UpperCount();  
  void Empty();
  void Clean();
  //  void DumpToFile(char FileName[100]);
  //  void ReadFromFile(char FileName[100]);
};

// void Point2Pattern::Print(){
//   long int i,j,k;
//   k = 0;
//   struct Point2 *TempCell;
//   for(i=0;i<=MaxXCell;i++){
//     for(j=0;j<=MaxYCell;j++){
//       //Rprintf("%d %d:\n",i,j);
//       TempCell = headCell[i][j]->next;
//       CHECK(TempCell, "internal error: TempCell is null in Print()");
// 	while(TempCell->next != TempCell){
// 	  k++;
// 	  Rprintf("%f %f %ld %ld %ld=%d %ld=%d UL0 %d UL1 %d %f\n",
// 		  TempCell->X,TempCell->Y,k,
// 		  TempCell->No,
// 		  i,int(TempCell->X/XCellDim),
// 		  j,int(TempCell->Y/YCellDim),
// 		  TempCell->InLower[0],TempCell->InLower[1],
// 		  TempCell->Beta);
// 	  TempCell = TempCell->next;
// 	  CHECK(TempCell, "internal error: TempCell is null in Print() loop");
//       }
//     }
//   }
//   Rprintf("Printed %ld points.\n",k);
// }

void Point2Pattern::Return(double *X, double *Y, int *num, int maxnum){
  long int i,j,k;
  k =0; *num = 0;
#ifdef DBGS
  Rprintf("executing Return()\n");
#endif
  if(UpperLiving[0]<=maxnum){
    struct Point2 *TempCell;
    for(i=0;i<=MaxXCell;i++){
      for(j=0;j<=MaxYCell;j++){
#ifdef DBGS
	//	Rprintf("%d %d:\n",i,j);
#endif
	TempCell = headCell[i][j]->next;
	CHECK(TempCell, "internal error: TempCell is null in Return()");
	while(TempCell->next != TempCell){
	  X[k] = TempCell->X;
	  Y[k] = TempCell->Y;	
	  k++;
	  TempCell = TempCell->next;
	  CHECK(TempCell, "internal error: TempCell is null in Return() loop");
	}
      }
    }    
    *num = k;
  } else {
    *num = -1;
  }
}

long int Point2Pattern::Count(){
  long int i,j,k;
  k = 0;
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      // Rprintf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      CHECK(TempCell, "internal error: TempCell is null in Count()");
      while(TempCell->next != TempCell){
	k++;
	TempCell = TempCell->next;
	CHECK(TempCell, "internal error: TempCell is null in Count() loop");
      }
    }
  }
  //Rprintf("Printed %d points.\n",k);
  return(k);
}

// a quick (over)estimate of the number of points in the pattern, 
// for storage allocation

long int Point2Pattern::UpperCount(){
  return(UpperLiving[0]);
}

void Point2Pattern::Empty(){
  struct Point2 *TempCell, *TempCell2;
  long int i,j;
  
#ifdef DBGS
  long int k;
  k=0;
  Rprintf("executing Empty()\n");
#endif

  for(i=0; i<=this->MaxXCell; i++){
    for(j=0; j<=this->MaxYCell; j++){
      TempCell = headCell[i][j]->next;
      CHECK(TempCell, "internal error: TempCell is null in Empty()");
      while(TempCell!=TempCell->next){	
#ifdef DBGS
	//	k++; Rprintf("%d %d %d\n",i,j,k);
#endif
	TempCell2 = TempCell->next;
	FREE(TempCell);
	TempCell = TempCell2;
	CHECK(TempCell, "internal error: TempCell is null in Empty() loop");
      }
      headCell[i][j]->next = dummyCell;
    }
  }
}

void Point2Pattern::Clean(){
  struct Point2 *TempCell, *TempCell2;
  long int i,j;
  
#ifdef DBGS
  Rprintf("executing Clean()\n");
#endif

  for(i=0; i<=MaxXCell; i++){
    for(j=0; j<=MaxYCell; j++){
      TempCell = headCell[i][j];
      CHECK(TempCell, "internal error: TempCell is null in Clean()");
      TempCell2 = headCell[i][j]->next;
      CHECK(TempCell2, "internal error: TempCell2 is null in Clean()");
      while(TempCell2!=TempCell2->next){
	TempCell2->No = 0;
	if(TempCell2->InLower[0]==0){
	  TempCell->next = TempCell2->next;
	  FREE(TempCell2);
	  TempCell2 = TempCell->next;
	  CHECK(TempCell2, 
		"internal error: TempCell2 is null in Clean() loop A");
	}
	else{
	  TempCell2 = TempCell2->next;
	  TempCell = TempCell->next;
	  CHECK(TempCell, "internal error: TempCell is null in Clean() loop B");
	  CHECK(TempCell2, 
		"internal error: TempCell2 is null in Clean() loop B");
	}
      }
    }
  }
}

//void Point2Pattern::DumpToFile(char FileName[100]){
//  FILE *out;
//  long int i,j;
//  out = fopen(FileName,"w");
//  struct Point2 *TempCell;
//  for(i=0;i<=MaxXCell;i++){
//    for(j=0;j<=MaxYCell;j++){
//    //Rprintf("%d %d:\n",i,j);
//    TempCell = headCell[i][j]->next;
//    while(TempCell->next != TempCell){
//	fprintf(out,"%f\t%f\t%ld\n",
//	       TempCell->X,TempCell->Y,TempCell->No);
//	TempCell = TempCell->next;
//    }
//  }
//}
//fclose(out);
//}

//void Point2Pattern::ReadFromFile(char FileName[100]){
//  FILE *out;
//long int k,XCell,YCell;
//float f1,xs,ys;
//out = fopen(FileName,"r");
//struct Point2 *TempCell;
//k=0;
//while(feof(out)==0){
//  k++;
//  fscanf(out,"%f%f\n",&xs,&ys);
//  //Rprintf("%f %f\n",xs,ys);
//  //
//  TempCell = ALLOCATE(struct Point2);
//  //
//  TempCell->No = k;
//  TempCell->X = xs;
//  TempCell->Y = ys;
//  TempCell->InLower[0] = 1;
//  TempCell->InLower[1] = 1;
//
//  f1 = (xs-Xmin)/XCellDim;  XCell = int(f1);
//  if(XCell>MaxXCell) XCell = MaxXCell;
//  f1 = (ys-Ymin)/YCellDim;  YCell = int(f1);
//  if(YCell>MaxYCell) YCell = MaxYCell;
//
//  TempCell->next = headCell[XCell][YCell]->next;
//  headCell[XCell][YCell]->next = TempCell;
//
//}
//fclose(out);
//Rprintf("%ld points loaded.\n",k);
//
//}


// ........................... Point processes ..........................
// ...................... (stationary, pairwise interaction) ............

class PointProcess {
 public:
  double Xmin, Xmax, Ymin, Ymax, TotalBirthRate, InteractionRange;
  PointProcess(double xmin, double xmax, double ymin, double ymax){
    Xmin = xmin; Xmax = xmax;
    Ymin = ymin; Ymax = ymax;
  }
  ~PointProcess(){}
  virtual void NewEvent(double *x, double *y, char *InWindow)=0;
  virtual void GeneratePoisson(Point *headPoint, 
			       long int *GeneratedPoints,
			       long int *LivingPoints,
			       long int *NoP)=0;
  virtual double Interaction(double dsquared)=0;
  //  virtual void CalcBeta(long int xsidepomm, long int ysidepomm, 
  //		   double *betapomm){ 
  //  Rprintf("Define CalcBeta...\n");
  // }
  //  virtual void CheckBeta(long int xsidepomm, long int ysidepomm, 
  //		   double *betapomm){ 
  //Rprintf("Define CheckBeta...\n");
  //}
  //  virtual double lnCondInt(struct Point2 *TempCell, Point2Pattern *p2p)
  //{ return(0.0);};
  //  virtual double lnDens(Point2Pattern *p2p);
  //  virtual void Beta(struct Point2 *TempCell){
  //    TempCell->Beta = 0;
  //    Rprintf("Define Beta...\n");};
};

//double PointProcess::lnDens(Point2Pattern *p2p){  
//// double f1;
//long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx;
//double dy,dx, lnDens,dst2;
//struct Point2 *TempCell, *TempCell2;
//
//dx = (Xmax-Xmin)/(double(p2p->MaxXCell+1));
//dy = (Ymax-Ymin)/(double(p2p->MaxYCell+1));
//rx = int(InteractionRange/dx+1.0);
//ry = int(InteractionRange/dy+1.0);
//
//  //Rprintf("1:%f 2:%f 3:%d 4:%d 5:%f 6:%f\n",dx,dy,rx,ry,
//  // this->InteractionRange,InteractionRange);
//  //Rprintf("mx:%d my:%d\n",p2p->MaxXCell,p2p->MaxYCell);
//
//  lnDens = 0;
//
//  //Rprintf("lnDens: %f (0)\n",lnDens);
//  
//  for(xc = 0; xc <= p2p->MaxXCell; xc++){
//    for(yc = 0; yc <= p2p->MaxYCell; yc++){
//      //if(xc==1) Rprintf("%d %d\n",xc,yc);
//      CHECK(p2p->headCell[xc][yc], 
//	    "internal error: p2p->headCell[xc][yc] is null in lnDens()");
//      TempCell = p2p->headCell[xc][yc]->next;
//      CHECK(TempCell, "internal error: TempCell is null in lnDens()");
//      while(TempCell != TempCell->next){
//	lnDens += log(TempCell->Beta);
//	//Rprintf("lnDens: %f (1) %d %d %d %d Beta %f\n",lnDens,xc,yc,
//	//       p2p->MaxXCell,p2p->MaxYCell,TempCell->Beta);
//	//if(lnDens<(-100000)){Rprintf("%f",lnDens); scanf("%f",&f1);}
//	if(InteractionRange>0){
//	  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
//	  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
//	  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
//	  if((yc-ry)>=0) fy=yc-ry; else fy = 0;
//	  for(xco = fx; xco <= tx; xco++){
//	    for(yco = fy; yco <= ty; yco++){
//	      //if(xc==1) Rprintf("%d %d %d %d %d %d\n",xco,yco,fx,tx,fy,ty);
//	      CHECK(p2p->headCell[xco][yco], 
//		    "internal error: p2p->headCell[xco][yco] is null in lnDens() loop");
//	      TempCell2 = p2p->headCell[xco][yco]->next;
//	      CHECK(TempCell2, 
//		    "internal error: TempCell2 is null in lnDens() loop A");
//	      while(TempCell2!=TempCell2->next){
//		if(TempCell2 != TempCell){
//		  dst2 = pow(TempCell->X-TempCell2->X,2)+
//			     pow(TempCell->Y-TempCell2->Y,2);
//		  lnDens += log(Interaction(dst2));
//		}
//		TempCell2 = TempCell2->next; 
//		CHECK(TempCell2, 
//		      "internal error: TempCell2 is null in lnDens() loop B");
//	      }
//	    }
//	  }
//	  //Rprintf("lnDens: %f\n",lnDens);
//	}
//	TempCell = TempCell->next;
//	CHECK(TempCell, 
//	      "internal error: TempCell is null in lnDens() at end");
//      }
//    }
//  }
//  return(lnDens);
//
//}

// ........................... Sampler ..........................

class Sampler{
 public:
  PointProcess *PP;
  Point2Pattern *P2P;
  long int GeneratedPoints, LivingPoints, NoP;
  //long int UpperLiving[2];
  Sampler(PointProcess *p){ PP = p;}
  ~Sampler(){}
  void Sim(Point2Pattern *p2p, long int *ST, long int *ET);
  long int BirthDeath(long int TimeStep,
		      struct Point *headLiving,
		      struct Point *headDeleted,
		      struct Point3 *headTransition);
  // WAS:  Sampler::Forward
  void Forward(long int TS, long int TT, char TX, char TY,
		      struct Point *Proposal, long int *DDD);
};


void Sampler::Forward(long int TS, long int TT, char TX, char TY,
		      struct Point *Proposal, long int *DDD){

  long int XCell, YCell, DirectionN;
  double dtmp2,dtmpx,dtmpy, tmpR, TempGamma[2], TempI;
  struct Point2 *TempCell, *TempCell2;
  float f1;

  /* Birth */
  if(TT==1){
    f1 = (Proposal->X-P2P->Xmin)/P2P->XCellDim;  XCell = int(f1);
    CLAMP(XCell, 0, P2P->MaxXCell, "XCell");
    f1 = (Proposal->Y-P2P->Ymin)/P2P->YCellDim;  YCell = int(f1);
    CLAMP(YCell, 0, P2P->MaxYCell, "YCell");
    //
    TempCell = ALLOCATE(struct Point2);
    //
    TempCell->No = Proposal->No;
    TempCell->X = Proposal->X;
    TempCell->Y = Proposal->Y;

    tmpR = Proposal->R;
    TempCell->next = P2P->headCell[XCell][YCell]->next;
    P2P->headCell[XCell][YCell]->next = TempCell;
    TempCell->InLower[0]=0;
    TempCell->InLower[1]=0;

    TempGamma[0] = 1.0; TempGamma[1] = 1.0;    

    /*same cell*/
    TempCell2 = TempCell->next; 
    CHECK(TempCell2, 
	  "internal error: TempCell2 is null in Forward() birth case");
    while(TempCell2 != TempCell2->next){
      dtmpx = TempCell->X - TempCell2->X;
      dtmpy = TempCell->Y - TempCell2->Y;
      dtmp2  = dtmpx*dtmpx+dtmpy*dtmpy;      
      TempI = PP->Interaction(dtmp2);
      if(TempCell2->InLower[0]==1) TempGamma[0] = TempGamma[0]*TempI;
      if(TempCell2->InLower[1]==1) TempGamma[1] = TempGamma[1]*TempI;
      TempCell2=TempCell2->next;
      CHECK(TempCell2, 
	    "internal error: TempCell2 is null in Forward() birth case loop");
    }
    /*eight other cells*/
    for(DirectionN=1;DirectionN<=8;DirectionN++){
      if(((XCell+P2P->DirX[DirectionN])>=0) &&
	 ((XCell+P2P->DirX[DirectionN])<=P2P->MaxXCell) &&
	 ((YCell+P2P->DirY[DirectionN])>=0) &&
	 ((YCell+P2P->DirY[DirectionN])<=P2P->MaxYCell)){
	CHECK(P2P->headCell[XCell+P2P->DirX[DirectionN]][YCell+P2P->DirY[DirectionN]], 
	      "internal error: HUGE P2P EXPRESSION is null in Forward() birth case loop A");
	TempCell2 = 
	  P2P->headCell[XCell+P2P->DirX[DirectionN]]
	  [YCell+P2P->DirY[DirectionN]]->next;
	CHECK(TempCell2, 
	      "internal error: TempCell2 is null in Forward() birth case loop B");
	while(TempCell2!=TempCell2->next){
	  dtmpx = TempCell->X - TempCell2->X;
	  dtmpy = TempCell->Y - TempCell2->Y;
	  dtmp2 = dtmpx*dtmpx+dtmpy*dtmpy;      
	  TempI = PP->Interaction(dtmp2);
	  if(TempCell2->InLower[0]==1) 
	    TempGamma[0] = TempGamma[0]*TempI;
	  if(TempCell2->InLower[1]==1) 
	    TempGamma[1] = TempGamma[1]*TempI;
	  TempCell2=TempCell2->next;
	CHECK(TempCell2, 
	      "internal error: TempCell2 is null in Forward() birth case loop C");
	}
      }
    }

    if(tmpR <= TempGamma[1] ){ 
      TempCell->InLower[0]=1;
      P2P->UpperLiving[0] = P2P->UpperLiving[0] +1;
    }
    if(tmpR <= TempGamma[0] ){ 
      TempCell->InLower[1]=1;
      P2P->UpperLiving[1] = P2P->UpperLiving[1] +1;
    }
  }
  /* Death */
  if(TT==0){
    TempCell=P2P->headCell[(int)TX][(int)TY];
    CHECK(TempCell, "internal error: TempCell is null in Forward() death case");
    while(TempCell->next->No != *DDD){
      TempCell = TempCell->next;
      CHECK(TempCell, 
	    "internal error: TempCell is null in Forward() death case loop");
      if(TempCell->next == TempCell) {
	// Rprintf("internal error: unexpected self-reference. Dumping...\n"); 
	// P2P->Print(); 
        error("internal error: unexpected self-reference");
	break;
      }
    };
    CHECK(TempCell->next, 
	  "internal error: TempCell->next is null in Forward() death case");
    if(*DDD!=TempCell->next->No) 
      Rprintf("diagnostic message: multi cell:  !!DDD:%ld TempUpper->No:%ld ",
	     *DDD,TempCell->No);
    if(TempCell->next->InLower[0]==1)
      P2P->UpperLiving[0] = P2P->UpperLiving[0] -1;
    if(TempCell->next->InLower[1]==1) 
      P2P->UpperLiving[1] = P2P->UpperLiving[1] -1;
    TempCell2 = TempCell->next;
    CHECK(TempCell2, 
	  "internal error: TempCell2 is null in Forward() death case B");
    TempCell->next = TempCell2->next;
    FREE(TempCell2);
    /* Common stuff */
    //KillCounter ++;
    *DDD = *DDD - 1;
  }
}


long int Sampler::BirthDeath(long int TimeStep,
		      struct Point *headLiving,
		      struct Point *headDeleted,
		      struct Point3 *headTransition){
  long int i,n;
  float f1,f2,f3,f4;
  double xtemp,ytemp;
  char InWindow, Success;
  struct Point *TempPoint, *TempPoint2;
  struct Point3 *TempTransition;

  R_CheckUserInterrupt();

  f1 = LivingPoints; f2 = PP->TotalBirthRate; f3 = f2/(f1+f2);
  f4 = slumptal();
  n = 0;
  Success = 0;

  //Rprintf("LivingPoints: %d TotalBirthRate %f GeneratedPoints %d\n",
  // LivingPoints,PP->TotalBirthRate,GeneratedPoints);
  
  /* Birth */
  while(Success==0){
  if(f4<f3){
    //Rprintf("Ping 1 (BD)\n");
    PP->NewEvent(&xtemp, &ytemp, &InWindow);
    //Rprintf("Ping 2 (BD)\n");
    if(InWindow==1){
      Success = 1;
      //
      TempTransition = ALLOCATE(struct Point3);
      //
      //Rprintf("Ping 3 (BD)\n");
      TempTransition->Case = 0;
      LivingPoints ++;
      GeneratedPoints ++;
      //
      TempPoint = ALLOCATE(struct Point);
      //
      TempPoint->X = xtemp;
      TempPoint->Y = ytemp;
      TempPoint->No = GeneratedPoints;
      TempPoint->R = slumptal();
      TempPoint->next = headLiving->next;
      headLiving->next = TempPoint;
      NoP ++;
      f1 = (TempPoint->X-P2P->Xmin)/P2P->XCellDim;  
      TempTransition->XCell = int(f1);
      f1 = (TempPoint->Y-P2P->Ymin)/P2P->YCellDim;  
      TempTransition->YCell = int(f1);
      
      //Rprintf("X %f XCell %d\n",TempPoint->X,TempTransition->XCell);
      // 
      CLAMP(TempTransition->XCell, 0, P2P->MaxXCell, "TempTransition->XCell");
      CLAMP(TempTransition->YCell, 0, P2P->MaxYCell, "TempTransition->YCell");
      TempTransition->next = headTransition->next;
      headTransition->next = TempTransition;
    }
  }
  /* Death */
  else{
    Success = 1;
    //
    TempTransition = ALLOCATE(struct Point3);
    //
    TempTransition->Case = 1;
    f1 = LivingPoints; f2 = f1*slumptal()+1.0;
    n = int(f2); if(n < 1) n = 1;
    if(n>LivingPoints){
      //      Rprintf("diagnostic message: random integer n=%ld > %ld = number of living points\n", n,LivingPoints);
      n=LivingPoints;
    }
    TempPoint2 = TempPoint = headLiving;
    for(i=1; i<=n; i++){ 
      TempPoint2 = TempPoint;
      TempPoint = TempPoint->next;
      }
    TempPoint2->next = TempPoint->next;
    
    TempPoint->next = headDeleted->next;  
    headDeleted->next = TempPoint;

    LivingPoints --;
    NoP --;
    TempTransition->next = headTransition->next;
    headTransition->next = TempTransition;
  }
  }
  return(n);
}

void Sampler::Sim(Point2Pattern *p2p, long int *ST, long int *ET) {

  P2P = p2p;
  long int StartTime, EndTime, TimeStep, D0Time, D0Living;
  long int XCell, YCell, DDD, i;
  float f1;
  
  /* Initialising linked listed for backward simulation */
  struct Point *headDeleted, *headLiving, *dummyDeleted, *dummyLiving;
  struct Point *TempPoint;
  //
  headLiving = ALLOCATE(struct Point);
  dummyLiving = ALLOCATE(struct Point);
  //
  headLiving->next = dummyLiving; dummyLiving->next = dummyLiving;
  //
  headDeleted = ALLOCATE(struct Point);
  dummyDeleted = ALLOCATE(struct Point);
  //
  headDeleted->next = dummyDeleted; dummyDeleted->next = dummyDeleted;

  struct Point2 *TempCell2;

  struct Point3 *headTransition, *dummyTransition;
  //
  headTransition = ALLOCATE(struct Point3);
  dummyTransition = ALLOCATE(struct Point3);
  //
  headTransition->next = dummyTransition; 
  dummyTransition->next = dummyTransition;
  
  PP->GeneratePoisson(headLiving, &GeneratedPoints,
			      &LivingPoints,
			      &NoP);  
    
  StartTime=1;
  EndTime=1;

  TimeStep = 0; D0Time = 0;
  D0Living = GeneratedPoints;

  long int tmp, D0;
  
  do{
    tmp=BirthDeath(TimeStep, headLiving,
		      headDeleted,
		      headTransition);
    if(tmp>0){ 
      if(tmp>(LivingPoints+1-D0Living)){
	D0Living --;
      }
    }
    D0Time++;
  }while(D0Living>0);
  tmp=BirthDeath(TimeStep, headLiving,
		      headDeleted,
		      headTransition); 
  StartTime=1; EndTime=D0Time+1; D0 = 0;

  do{	 
    if(D0==1){
      for(TimeStep=StartTime;TimeStep<=EndTime;TimeStep ++){
	tmp=BirthDeath(TimeStep, headLiving,
		       headDeleted,
		       headTransition);      
      }
    }
    D0 = 1;
    P2P->Empty();
    
    /*
    headUpper->next = dummyUpper; dummyUpper->next = dummyUpper;
    for(XCell=0;XCell<=P2P->MaxXCell;XCell++){
      for(YCell=0;YCell<=P2P->MaxYCell;YCell++){
	headUpperCell[XCell][YCell]->next=dummyUpper;
      }
    }
    */
    
    P2P->UpperLiving[0] = LivingPoints;
    P2P->UpperLiving[1] = 0;
    
    P2P->NoP = 0;
    i=0;
    TempPoint = headLiving->next;
    CHECK(TempPoint, "internal error: TempPoint is null in Sim()");
    while(TempPoint!=TempPoint->next){
      i++;
      //
      TempCell2 = ALLOCATE(struct Point2);
      //
      TempCell2->No = TempPoint->No;
      TempCell2->X = TempPoint->X;
      TempCell2->Y = TempPoint->Y;
      TempCell2->InLower[0] = 1;
      TempCell2->InLower[1] = 0;
      f1 = (TempPoint->X-P2P->Xmin)/P2P->XCellDim;  XCell = int(floor(f1));
      CLAMP(XCell, 0, P2P->MaxXCell, "XCell");
      f1 = (TempPoint->Y-P2P->Ymin)/P2P->YCellDim;  YCell = int(floor(f1));
      CLAMP(YCell, 0, P2P->MaxYCell, "YCell");
      TempCell2->next = P2P->headCell[XCell][YCell]->next;
      P2P->headCell[XCell][YCell]->next = TempCell2;
      
      TempPoint = TempPoint->next;
      CHECK(TempPoint, "internal error: TempPoint is null in Sim() loop");
    }
    
    //P2P->DumpToFile("temp0.dat");
    
    struct Point3 *TempTransition;
    struct Point *Proposal;
    
    TempTransition = headTransition->next;
    CHECK(TempTransition, "internal error: TempTransition is null in Sim()");
    Proposal = headDeleted->next;
    DDD = GeneratedPoints;
    
    for(TimeStep=EndTime;TimeStep>=1;TimeStep--){
      R_CheckUserInterrupt();
      Forward(TimeStep,TempTransition->Case,
	      TempTransition->XCell,TempTransition->YCell,
	      Proposal,&DDD);
      if(TempTransition->Case == 1) Proposal = Proposal ->next;
      TempTransition = TempTransition->next;
      CHECK(TempTransition, 
	    "internal error: TempTransition is null in Sim() loop");
    }
    
    /* Doubling strategy used!*/
    StartTime = EndTime+1;
    EndTime=EndTime*2;
    
    //P2P->DumpToFile("temp.dat");
    
  }while(P2P->UpperLiving[0]!=P2P->UpperLiving[1]);
  P2P->Clean();
  i=0;
  struct Point *TempPoint2;
  TempPoint = headLiving;
  TempPoint2 = headLiving->next;
  CHECK(TempPoint2, 
	"internal error: TempPoint2 is null in Sim() position B");
  while(TempPoint!=TempPoint->next){
    i++;
    FREE(TempPoint);
    TempPoint = TempPoint2;
    TempPoint2 = TempPoint2->next;
    CHECK(TempPoint2, 
	  "internal error: TempPoint2 is null in Sim() loop C");
  }
  FREE(TempPoint);
  
  i = 0;
  TempPoint = headDeleted;
  TempPoint2 = headDeleted->next;
  CHECK(TempPoint2, 
	"internal error: TempPoint2 is null in Sim() position D");
  while(TempPoint!=TempPoint->next){
    i++;
    FREE(TempPoint);
    TempPoint = TempPoint2;
    TempPoint2 = TempPoint2->next;
    CHECK(TempPoint2, 
	  "internal error: TempPoint2 is null in Sim() loop D");
  }
  FREE(TempPoint);
  //Rprintf("%d ",i);

  struct Point3 *TempTransition,*TempTransition2;

  i = 0;
  TempTransition = headTransition;
  TempTransition2 = headTransition->next;
  CHECK(TempTransition2, 
	"internal error: TempTransition2 is null in Sim() position E");
  while(TempTransition!=TempTransition->next){
    i++;
    FREE(TempTransition);
    TempTransition = TempTransition2;
    TempTransition2 = TempTransition2->next;
    CHECK(TempTransition2, 
	  "internal error: TempTransition2 is null in Sim() loop F");
  }
  FREE(TempTransition);
  //Rprintf("%d ST: %d ET: %d\n",i,StartTime,EndTime);
  //scanf("%f",&f1);
  *ST = StartTime;
  *ET = EndTime;
}

#include "PerfectStrauss.h"
#include "PerfectStraussHard.h"
#include "PerfectHardcore.h"
#include "PerfectDiggleGratton.h"
#include "PerfectDGS.h"
#include "PerfectPenttinen.h"
