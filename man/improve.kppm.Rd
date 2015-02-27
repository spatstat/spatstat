\name{improve.kppm}
\alias{improve.kppm}
\title{Improve Intensity Estimate of Fitted Cluster Point Process Model}
\description{
  Update the fitted intensity of a fitted cluster point process model.
}
\usage{
improve.kppm(object, type=c("quasi", "wclik1", "clik1"), rmax = NULL,
             eps.rmax = 0.01, dimyx = 50, maxIter = 100, tolerance = 1e-06,
             fast = TRUE, vcov = FALSE, fast.vcov = FALSE, verbose = FALSE,
                          save.internals = FALSE)
}
\arguments{
  \item{object}{
    Fitted cluster point process model (object of class \code{"kppm"}).
  }
  \item{type}{
    A character string indicating the method of estimation.
    Current options are \code{"clik1"}, \code{"wclik1"} and \code{"quasi"}
    for, respectively, first order composite (Poisson) likelihood,
    weighted first order composite likelihood and quasi-likelihood.
  }
  \item{rmax}{
    Optional. The dependence range. Not usually specified by the user.
  }
  \item{eps.rmax}{
    Numeric. A small positive number which is used to determine \code{rmax}
    from the tail behaviour of the pair correlation function. Namely
    \code{rmax} is the smallest value of \eqn{r}
    at which \eqn{(g(r)-1)/(g(0)-1)}
    falls below \code{eps.rmax}. 
    Ignored if \code{rmax} is provided.
  }
  \item{dimyx}{
    Pixel array dimensions. See Details.
  }
  \item{maxIter}{
    Integer. Maximum number of iterations of iterative weighted least squares
    (Fisher scoring).
  }
  \item{tolerance}{
    Numeric. Tolerance value specifying when to stop iterative weighted
    least squares (Fisher scoring).
  }
  \item{fast}{
    Logical value indicating whether tapering should be used to make the 
    computations faster (requires the package \pkg{Matrix}).
  }
  \item{vcov}{
    Logical value indicating whether to calculate the asymptotic variance
    covariance/matrix.
  }
  \item{fast.vcov}{
    Logical value indicating whether tapering should be used for the
    variance/covariance matrix to make the computations faster
    (requires the package \pkg{Matrix}). Caution:
    This is expected to underestimate the true asymptotic variances/covariances.
  }
  \item{verbose}{
    A logical indicating whether the details of computations should be printed.
  }
  \item{save.internals}{
    A logical indicating whether internal quantities should be saved in the
    returned object (mostly for development purposes).
  }
}
\value{
  A fitted cluster point process model of class \code{"kppm"}.
}
\details{
  This function reestimates the intensity parameters in a fitted \code{"kppm"}
  object. If \code{type="clik1"} estimates are based on the first order
  composite (Poisson) likelihood, which ignores dependence between the
  points. Note that \code{type="clik1"} is mainly included for testing
  purposes and is not recommended for the typical user;
  instead the more efficient \code{\link{kppm}}
  with \code{improve.type="none"} should be used.
  
  When \code{type="quasi"} or \code{type="wclik1"} the dependence
  structure between the points is incorporated in the estimation
  procedure by using the estimated pair correlation function in the
  estimating equation.

  In all cases the estimating equation is based on dividing the
  observation window into small subregions and count the number of points
  in each subregion. To do this the observation window is first
  converted into a digital mask by \code{\link{as.mask}} where the
  resolution is controlled by the argument \code{dimyx}. The
  computational time grows with the cube of the number of subregions, so fine
  grids may take very long to compute (or even run out of memory).
}
\seealso{
  \code{\link{ppm}},
  \code{\link{kppm}},
  \code{\link{improve.kppm}}
}
\references{
  Waagepetersen, R. (2007) An estimating function approach to inference
  for inhomogeneous Neyman-Scott processes, \emph{Biometrics},
  \bold{63}, 252-258.
  
  Guan, Y. and Shen, Y. (2010) A weighted estimating equation approach
  to inference for inhomogeneous spatial point processes, \emph{Biometrika},
  \bold{97}, 867-880.
  
  Guan, Y., Jalilian, A. and Waagepetersen, R. (2013) 
  Quasi-likelihood for spatial Cox point processes. 
  \emph{Submitted}.
}
\examples{
  # fit a Thomas process using minimum contrast estimation method 
  # to model interaction between points of the pattern
  fit0 <- kppm(bei ~ elev + grad, data = bei.extra)

  # fit the log-linear intensity model with quasi-likelihood method
  fit1 <- improve.kppm(fit0, type="quasi")

  # compare
  coef(fit0)
  coef(fit1)
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
%  \url{http://www.razi.ac.ir/ajalilian/}
  and Rasmus Waagepetersen
  \email{rw@math.aau.dk}
  adapted for \pkg{spatstat} by Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{fit model}

