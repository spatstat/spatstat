\name{exactMPLEstrauss}
\alias{exactMPLEstrauss}
\title{
  Exact Maximum Pseudolikelihood Estimate for Stationary Strauss Process
}
\description{
  Computes, to very high accuracy, the Maximum Pseudolikelihood Estimates
  of the parameters of a stationary Strauss point process.
}
\usage{
  exactMPLEstrauss(X, R, ngrid = 2048, plotit = FALSE, project=TRUE)
}
\arguments{
  \item{X}{
    Data to which the Strauss process will be fitted.
    A point pattern dataset (object of class \code{"ppp"}).
  }
  \item{R}{
    Interaction radius of the Strauss process.
    A non-negative number.
  }
  \item{ngrid}{
    Grid size for calculation of integrals. An integer, giving the
    number of grid points in the \eqn{x} and \eqn{y} directions.
  }
  \item{plotit}{
    Logical. If \code{TRUE}, the log pseudolikelihood is plotted
    on the current device. 
  }
  \item{project}{
    Logical. If \code{TRUE} (the default), the parameter
    \eqn{\gamma}{gamma} is constrained to lie in the interval
    \eqn{[0,1]}. If \code{FALSE}, this constraint is not applied.
  }
}
\details{
  This function is intended mainly for technical investigation
  of algorithm performance. Its practical use is quite limited.

  It fits the stationary Strauss point process model
  to the point pattern dataset \code{X} by maximum pseudolikelihood
  (with the border edge correction) using an algorithm with very high accuracy.
  This algorithm is more accurate than the
  \emph{default} behaviour of the model-fitting function
  \code{\link{ppm}} because the discretisation is much finer.

  Ripley (1988) and Baddeley and Turner (2000) derived the
  log pseudolikelihood for the stationary Strauss
  process, and eliminated the parameter \eqn{\beta}{beta},
  obtaining an exact formula for the partial log pseudolikelihood
  as a function of the interaction parameter \eqn{\gamma}{gamma} only.
  The algorithm evaluates this expression to a high degree of accuracy,
  using numerical integration on a \code{ngrid * ngrid} lattice,
  uses \code{\link[stats]{optim}} to maximise the log pseudolikelihood
  with respect to \eqn{\gamma}{gamma}, and finally recovers
  \eqn{\beta}{beta}.

  The result is a vector of length 2, containing the fitted coefficients
  \eqn{\log\beta}{log(beta)} and \eqn{\log\gamma}{log(gamma)}.
  These values correspond to the entries that would be obtained with
  \code{coef(ppm(X, ~1, Strauss(R)))}.
  The fitted coefficients are typically accurate to
  within \eqn{10^{-6}}{10^(-6)} as shown in Baddeley and Turner (2013).
  
  Note however that (by default) \code{exactMPLEstrauss} 
  constrains the parameter \eqn{\gamma}{gamma} to lie in the
  interval \eqn{[0,1]} in which the point process is well defined
  (Kelly and Ripley, 1976)
  whereas \code{\link{ppm}} does not constrain
  the value of \eqn{\gamma}{gamma} (by default). This behaviour is controlled by
  the argument \code{project} to \code{\link{ppm}} and
  \code{exactMPLEstrauss}. The default for \code{\link{ppm}}
  is \code{project=FALSE}, while the default for \code{exactMPLEstrauss}
  is \code{project=TRUE}.
}
\value{
  Vector of length 2.
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.

  Baddeley, A. and Turner, R. (2013)
  Manuscript in preparation.  

  Kelly, F.P. and Ripley, B.D. (1976)
  On Strauss's model for clustering.
  \emph{Biometrika} \bold{63}, 357--360.

  Ripley, B.D. (1988)
  \emph{Statistical inference for spatial processes}.
  Cambridge University Press.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{ppm}}
}
\examples{
\testonly{
   exactMPLEstrauss(cells, 0.1, ngrid=128)
   exactMPLEstrauss(cells, 0.1, ngrid=128, project=FALSE)
}
if(interactive()) {
   exactMPLEstrauss(cells, 0.1)
   coef(ppm(cells, ~1, Strauss(0.1)))
   coef(ppm(cells, ~1, Strauss(0.1), nd=128))
   exactMPLEstrauss(redwood, 0.04)
   exactMPLEstrauss(redwood, 0.04, project=FALSE)
   coef(ppm(redwood, ~1, Strauss(0.04)))
}
}
\keyword{spatial}
\keyword{models}
