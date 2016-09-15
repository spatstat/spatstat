\name{bc.ppm}
\alias{bc}
\alias{bc.ppm}
\title{
  Bias Correction for Fitted Model
}
\description{
  Applies a first-order bias correction to a fitted model.
}
\usage{
  bc(fit, \dots)

  \method{bc}{ppm}(fit, \dots, nfine = 256)
}
\arguments{
  \item{fit}{
    A fitted point process model (object of class \code{"ppm"})
    or a model of some other class.
  }
  \item{\dots}{
    Additional arguments are currently ignored.
  }
  \item{nfine}{
    Grid dimensions for fine grid of locations.
    An integer, or a pair of integers. See Details.
  }
}
\details{
  This command applies the first order Newton-Raphson bias correction method of
  Baddeley and Turner (2014, sec 4.2) to a fitted model.
  The function \code{bc} is generic, with a method for fitted point
  process models of class \code{"ppm"}.

  A fine grid of locations, of dimensions \code{nfine * nfine} or
  \code{nfine[2] * nfine[1]}, is created over the original window of the
  data, and the intensity or conditional intensity of the fitted model is
  calculated on this grid. The result is used to update the fitted
  model parameters once by a Newton-Raphson update.

  This is only useful if the quadrature points used to fit the original
  model \code{fit} are coarser than the grid of points specified by
  \code{nfine}. 
}
\value{
  A numeric vector, of the same length as \code{coef(fit)}, giving updated
  values for the fitted model coefficients.
}
\references{
  Baddeley, A. and Turner, R. (2014)
  Bias correction for parameter estimates of spatial point process models.
  \emph{Journal of Statistical Computation and Simulation}
  \bold{84}, 1621--1643.
  DOI: 10.1080/00949655.2012.755976
}
\author{
  \adrian
  and
  \rolf.
}
\seealso{
  \code{\link{rex}}
}
\examples{
  fit <- ppm(cells ~ x, Strauss(0.07))
  coef(fit)
  if(!interactive()) {
    bc(fit, nfine=64)
  } else {
    bc(fit)
  }
}
\keyword{spatial}
\keyword{models}
\keyword{math}
\keyword{optimize}
