\name{rex}
\alias{rex}
\title{
  Richardson Extrapolation
}
\description{
  Performs Richardson Extrapolation on a sequence of
  approximate values.
}
\usage{
rex(x, r = 2, k = 1, recursive = FALSE)
}
\arguments{
  \item{x}{
    A numeric vector or matrix,
    whose columns are successive estimates or approximations
    to a vector of parameters.
  }
  \item{r}{
    A number greater than 1.
    The ratio of successive step sizes.
    See Details.
  }
  \item{k}{
    Integer. The order of convergence assumed. See Details.
  }
  \item{recursive}{
    Logical value indicating whether to perform one step of
    Richardson extrapolation (\code{recursive=FALSE}, the default)
    or repeat the extrapolation procedure until a best estimate
    is obtained (\code{recursive=TRUE}.
  }
}
\details{
  Richardson extrapolation is a general technique for improving numerical
  approximations, often used in numerical integration (Brezinski and
  Zaglia, 1991). It can also be used to improve parameter estimates
  in statistical models (Baddeley and Turner, 2014). 
  
  The successive columns of \code{x} are assumed to
  have been obtained using approximations with step sizes
  \eqn{a, a/r, a/r^2, \ldots}{a, a/r, a/r^2, ...}
  where \eqn{a} is the initial step size (which does not need to be
  specified).

  Estimates based on a step size \eqn{s} are assumed to have an error
  of order \eqn{s^k}.

  Thus, the default values \code{r=2} and {k=1} imply that the errors in
  the second column of \code{x} should be roughly \eqn{(1/r)^k = 1/2} as large
  as the errors in the first column, and so on.
}
\value{
  A matrix whose columns contain a sequence of improved estimates.
}
\references{
  Baddeley, A. and Turner, R. (2014)
  Bias correction for parameter estimates of spatial point process models.
  \emph{Journal of Statistical Computation and Simulation}
  \bold{84}, 1621--1643.
  DOI: 10.1080/00949655.2012.755976

  Brezinski, C. and Zaglia, M.R. (1991)
  \emph{Extrapolation Methods. Theory and Practice}.
  North-Holland.
}
\author{
  \adrian
  and
  \rolf.
}
\seealso{
  \code{\link{bc}}
}
\examples{
   # integrals of sin(x) and cos(x) from 0 to pi
   # correct answers: 2, 0
   est <- function(nsteps) {
     xx <- seq(0, pi, length=nsteps)
     ans <- pi * c(mean(sin(xx)), mean(cos(xx)))
     names(ans) <- c("sin", "cos")
     ans
   }
   X <- cbind(est(10), est(20), est(40))
   X
   rex(X)
   rex(X, recursive=TRUE)

   # fitted Gibbs point process model
   fit0 <- ppm(cells ~ 1, Strauss(0.07), nd=16)
   fit1 <- update(fit0, nd=32)
   fit2 <- update(fit0, nd=64)
   co <- cbind(coef(fit0), coef(fit1), coef(fit2))
   co 
   rex(co, k=2, recursive=TRUE)
}
\keyword{math}
\keyword{optimize}

