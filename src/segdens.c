#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

/*
  segdens.c

  Convolution of segments with Gaussian kernel

  Adrian Baddeley, 02 dec 2016 
  Licence: GPL >= 2.0
*/


#define DNORM(X, SIG) dnorm((X), (double) 0.0, (SIG), FALSE)

#define PNORM(X, SIG) pnorm((X), (double) 0.0, (SIG), TRUE, FALSE)

void segdens(sigma, ns, xs, ys, alps, lens, np, xp, yp, z) 
     double *sigma; /* bandwidth */
     int *ns; /* number of line segments */
     double *xs, *ys, *alps, *lens;  /* first endpoint, angle, length */
     int *np; /* number of pixels or test locations */
     double *xp, *yp; /* pixel coordinates */
     double *z; /* result, assumed initially 0 */
{
  int i, j, Ns, Np;
  double Sigma;
  double xsi, ysi, angi, leni, cosi, sini;
  double dx, dy, u1, u2;

  Ns = *ns;
  Np = *np;
  Sigma = *sigma;

  for(i = 0; i < Ns; i++) {
    R_CheckUserInterrupt();
    xsi = xs[i];
    ysi = ys[i];
    angi = alps[i];
    leni = lens[i];
    cosi = cos(angi);
    sini = sin(angi);
    for(j = 0; j < Np; j++) {
      dx = xp[j] - xsi;
      dy = yp[j] - ysi;
      u1 = dx * cosi + dy * sini;
      u2 = -dx * sini + dy * cosi;
      z[j] += DNORM(u2, Sigma) * (PNORM(u1, Sigma) - PNORM(u1-leni, Sigma));
    }
  }
}
