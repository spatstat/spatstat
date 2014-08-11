#include <Rmath.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"
#include "crossloop.h"
#include "constants.h"
/*

  densptcross.c

  $Revision: 1.2 $     $Date: 2014/04/02 10:27:43 $

  Assumes point patterns are sorted in increasing order of x coordinate

  *crdenspt     Density estimate at points
  *crsmoopt     Smoothed mark values at points

*/

#define TWOPI M_2PI

double sqrt(), exp();

#define STD_DECLARATIONS				\
  int i, j, n1, n2, maxchunk, jleft;                    \
  double x1i, y1i, xleft, dx, dy, d2, rmax, r2max;      \
  double *x1, *y1, *x2, *y2;

#define STD_INITIALISE				\
  n1 = *nquery;					\
  x1 = xq; y1 = yq;                             \
  n2 = *ndata;					\
  x2 = xd; y2 = yd;                             \
  rmax = *rmaxi;				\
  r2max = rmax * rmax


/* ----------------- density estimation -------------------- */

void crdenspt(nquery, xq, yq, ndata, xd, yd, rmaxi, sig, result) 
     /* inputs */
     int *nquery;            /* number of locations to be interrogated */
     double *xq, *yq;    /* (x,y) coordinates to be interrogated */
     int *ndata;            /* number of data points */
     double *xd, *yd;    /* (x,y) coordinates of data */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *sig;      /* Gaussian sd */
     /* output */
     double *result;   /* vector of computed density values */
{
  STD_DECLARATIONS;
  double resulti, coef;	
  double sigma, twosig2; 
  STD_INITIALISE;

  sigma = *sig;				      
  twosig2 = 2.0 * sigma * sigma;	
  coef = 1.0/(TWOPI * sigma * sigma);  

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP( { resulti = 0.0; },
            { resulti += exp(-d2/twosig2); } ,
	    { result[i] = coef * resulti; })

}


void wtcrdenspt(nquery, xq, yq, ndata, xd, yd, wd, rmaxi, sig, result) 
     /* inputs */
     int *nquery;        /* number of locations to be interrogated */
     double *xq, *yq;    /* (x,y) coordinates to be interrogated */
     int *ndata;         /* number of data points */
     double *xd, *yd;    /* (x,y) coordinates of data */
     double *wd;         /* weights of data points */
     double *rmaxi;      /* maximum distance at which points contribute */
     double *sig;        /* Gaussian sd */
     /* output */
     double *result;   /* vector of computed density values */
{
  STD_DECLARATIONS;
  double resulti, coef;	
  double sigma, twosig2; 
  STD_INITIALISE;

  sigma = *sig;				      
  twosig2 = 2.0 * sigma * sigma;	
  coef = 1.0/(TWOPI * sigma * sigma);  

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP( { resulti = 0.0; },
	    { resulti += wd[j] * exp(-d2/twosig2); },
	    { result[i] = coef * resulti; } )

 }

/* ------------- anisotropic versions -------------------- */

void acrdenspt(nquery, xq, yq, ndata, xd, yd, rmaxi, detsigma, sinv, result) 
     /* inputs */
     int *nquery;            /* number of locations to be interrogated */
     double *xq, *yq;    /* (x,y) coordinates to be interrogated */
     int *ndata;            /* number of data points */
     double *xd, *yd;    /* (x,y) coordinates of data */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *detsigma;  /* determinant of variance matrix */
     double *sinv;      /* inverse variance matrix (2x2, flattened) */
     /* output */
     double *result;   /* vector of computed density values */
{
  STD_DECLARATIONS;
  double resulti, coef;	
  double detsig, s11, s12, s21, s22;
  STD_INITIALISE;
  detsig = *detsigma;
  coef = 1.0/(TWOPI * sqrt(detsig));
  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP( { resulti = 0.0; },
	    { resulti += exp(-(dx * (dx * s11 + dy * s12) \
			       + dy * (dx * s21 + dy * s22))/2.0); },
	    { result[i] = coef * resulti; })
}


void awtcrdenspt(nquery, xq, yq, ndata, xd, yd, wd, rmaxi, detsigma, sinv, result) 
     /* inputs */
     int *nquery;        /* number of locations to be interrogated */
     double *xq, *yq;    /* (x,y) coordinates to be interrogated */
     int *ndata;         /* number of data points */
     double *xd, *yd;    /* (x,y) coordinates of data */
     double *wd;         /* weights of data points */
     double *rmaxi;      /* maximum distance at which points contribute */
     double *detsigma;   /* determinant of variance matrix */
     double *sinv;       /* inverse variance matrix (2x2, flattened) */
     /* output */
     double *result;     /* vector of weighted density values */
{
  STD_DECLARATIONS;
  double resulti, coef;	
  double detsig, s11, s12, s21, s22;
  STD_INITIALISE;
  detsig = *detsigma;
  coef = 1.0/(TWOPI * sqrt(detsig));
  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP( { resulti = 0.0; },
	    { resulti += wd[j] * \
		exp(-(dx * (dx * s11 + dy * s12)			\
		      + dy * (dx * s21 + dy * s22))/2.0); },
	    { result[i] = coef * resulti; })
 }


/* --------------- smoothing --------------------------- */

void crsmoopt(nquery, xq, yq, ndata, xd, yd, vd, rmaxi, sig, result) 
     /* inputs */
     int *nquery;            /* number of locations to be interrogated */
     double *xq, *yq;    /* (x,y) coordinates to be interrogated */
     int *ndata;            /* number of data points */
     double *xd, *yd;    /* (x,y) coordinates of data */
     double *vd;         /* mark values at data points */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *sig;      /* Gaussian sd */
     /* output */
     double *result;   /* vector of computed smoothed values */
{
  STD_DECLARATIONS;
  double sigma, twosig2;
  double numer, denom, wij; 

  STD_INITIALISE;
  sigma = *sig;
  twosig2 = 2.0 * sigma * sigma;

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP({ numer = denom = 0.0; },
	   { \
	     wij = exp(-d2/twosig2);		\
	     denom += wij;			\
	     numer += wij * vd[j];		\
	   },					
	   {					\
	     result[i] = numer/denom;		\
	   })
 }


void wtcrsmoopt(nquery, xq, yq, ndata, xd, yd, vd, wd, rmaxi, sig, result) 
     /* inputs */
     int *nquery;            /* number of locations to be interrogated */
     double *xq, *yq;    /* (x,y) coordinates to be interrogated */
     int *ndata;            /* number of data points */
     double *xd, *yd;    /* (x,y) coordinates of data */
     double *vd;         /* mark values at data points */
     double *wd;         /* weights of data points */
     double *rmaxi;    /* maximum distance */
     double *sig;      /* Gaussian sd */
     /* output */
     double *result;    /* vector of computed smoothed values */
{
  STD_DECLARATIONS;
  double sigma, twosig2;
  double numer, denom, wij; 

  STD_INITIALISE;
  sigma = *sig;
  twosig2 = 2.0 * sigma * sigma;

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP({ numer = denom = 0.0; },
	   {						\
	     wij = wd[j] * exp(-d2/twosig2);	\
	     denom += wij;				\
	     numer += wij * vd[j];			\
	   },						
	   {						\
	     result[i] = numer/denom;			\
	   })
}

/* ------------- anisotropic versions -------------------- */

void acrsmoopt(nquery, xq, yq, ndata, xd, yd, vd, rmaxi, sinv, result) 
     /* inputs */
     int *nquery;            /* number of locations to be interrogated */
     double *xq, *yq;    /* (x,y) coordinates to be interrogated */
     int *ndata;            /* number of data points */
     double *xd, *yd;    /* (x,y) coordinates of data */
     double *vd;         /* mark values at data points */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *sinv;      /* inverse variance matrix (2x2, flattened) */
     /* output */
     double *result;   /* vector of smoothed values */
{
  STD_DECLARATIONS;
  double s11, s12, s21, s22;
  double numer, denom, wij; 

  STD_INITIALISE;
  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP({ numer = denom = 0.0; },
	   {							\
	     wij = exp(-(dx * (dx * s11 + dy * s12)		\
			 + dy * (dx * s21 + dy * s22))/2.0);	\
	     denom += wij;					\
	     numer += wij * vd[j];				\
	   },
	   {					\
	     result[i] = numer/denom;		\
	   })
}


void awtcrsmoopt(nquery, xq, yq, ndata, xd, yd, vd, wd, rmaxi, sinv, result) 
     /* inputs */
     int *nquery;        /* number of locations to be interrogated */
     double *xq, *yq;    /* (x,y) coordinates to be interrogated */
     int *ndata;         /* number of data points */
     double *xd, *yd;    /* (x,y) coordinates of data */
     double *vd;         /* mark values at data points */
     double *wd;         /* weights of data points */
     double *rmaxi;      /* maximum distance at which points contribute */
     double *sinv;       /* inverse variance matrix (2x2, flattened) */
     /* output */
     double *result;    /* vector of smoothed values */
{
  STD_DECLARATIONS;
  double s11, s12, s21, s22;
  double numer, denom, wij; 

  STD_INITIALISE;

  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP({ numer = denom = 0.0; },
	   {								\
	     wij = wd[j] * exp(-(dx * (dx * s11 + dy * s12)		\
				     + dy * (dx * s21 + dy * s22))/2.0); \
	     denom += wij;						\
	     numer += wij * vd[j];					\
	   },
	   {					\
	     result[i] = numer/denom;		\
	   })
}

