/* 
   dist2.h 

   External declarations for the functions defined in dist2.c
   and
   In-line cpp macros for similar purposes

   $Revision: 1.20 $ $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

double dist2(double u, double v, double x, double y, double *period);

double dist2either(double u, double v, double x, double y, double *period);

int dist2thresh(double u, double v, double x, double y, double *period, double r2);

int dist2Mthresh(double u, double v, double x, double y, double *period, double r2);

/* 
   Efficient macros to test closeness of points
*/

/* 
   These variables must be declared
   (note: some files e.g. straush.c use 'RESIDUE' explicitly)
*/

#define DECLARE_CLOSE_VARS \
  register double DX, DY, DXP, DYP, RESIDUE

#define DECLARE_CLOSE_D2_VARS \
  register double DX, DY, DXP, DYP, DX2

#define CLOSE(U,V,X,Y,R2)		\
  ((DX = X - U),			\
   (RESIDUE = R2 - DX * DX),		\
   ((RESIDUE > 0.0) &&			\
    ((DY = Y - V),                      \
     (RESIDUE = RESIDUE - DY * DY),     \
     (RESIDUE > 0.0))))

#define CLOSE_D2(U,V,X,Y,R2,D2)						\
  ((DX = X - U),							\
   (DX2 = DX * DX),							\
   (DX2 < R2) && (((DY = Y - V),					\
		   (D2 = DX2 + DY * DY),				\
		   (D2 < R2))))

/*
  The following calculates X mod P, 
  but it works only if X \in [-P, P]
  so that X is the difference between two values
  that lie in an interval of length P 
*/

#define CLOSE_PERIODIC(U,V,X,Y,PERIOD,R2)				\
  ((DX = X - U),							\
   (DX = (DX < 0.0) ? -DX : DX),					\
   (DXP = (PERIOD)[0] - DX),						\
   (DX = (DX < DXP) ? DX : DXP),					\
   (RESIDUE = R2 - DX * DX),						\
   ((RESIDUE > 0.0) && ((DY = Y - V),					\
			(DY = (DY < 0.0) ? -DY : DY),			\
			(DYP = (PERIOD)[1] - DY),			\
			(DY = (DY < DYP) ? DY : DYP),			\
                        (RESIDUE = RESIDUE - DY * DY),                  \
			(RESIDUE > 0.0) )))

#define CLOSE_PERIODIC_D2(U,V,X,Y,PERIOD,R2,D2)				\
  ((DX = X - U),							\
   (DX = (DX < 0.0) ? -DX : DX),					\
   (DXP = (PERIOD)[0] - DX),						\
   (DX = (DX < DXP) ? DX : DXP),					\
   (D2 = DX * DX),							\
   ((D2 < R2) && ((DY = Y - V),						\
		  (DY = (DY < 0.0) ? -DY : DY),				\
		  (DYP = (PERIOD)[1] - DY),				\
		  (DY = (DY < DYP) ? DY : DYP),				\
		  (D2 += DY * DY),					\
		  (D2 < R2) )))





