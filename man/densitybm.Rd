\name{densitybm}
\alias{densitybm}
\title{
  Barry-McIntyre Density Estimate of Point Pattern
}
\description{
  Computes the Barry-McIntyre estimate of the density of a point pattern.
}
\usage{
densitybm(X, sigma, eps = NULL)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{sigma}{
    Smoothing bandwidth. A single number giving the equivalent
    standard deviation of the smoother.
  }
  \item{eps}{
    Optional. A single number specifying the width and height of pixels.
  }
}
\details{
  Barry and McIntyre (2011) proposed an alternative to the usual
  kernel density (or intensity) estimator for point pattern data observed in an
  irregular spatial window. Their estimator is based on a random walk
  on the pixel grid inside the observation window. It is automatically
  edge-corrected.
}
\value{
  Pixel image (object of class \code{"im"}) giving the estimated
  intensity of the point process.
}
\seealso{
  \code{\link{density.ppp}} for the usual kernel estimator,
  and \code{\link{adaptive.density}} for the
  tessellation-based estimator.
}
\references{
  Barry, R.P. and McIntyre, J. (2011)
  Estimating animal densities and home range in regions with irregular
  boundaries and holes: a lattice-based alternative to the kernel
  density estimator. \emph{Ecological Modelling} \bold{222}, 1666--1672.
}
\author{
  \adrian
  and Greg McSwiggan.
}
\examples{
   X <- runifpoint(25, letterR)
   Z <- densitybm(X, 0.2)
   plot(Z, main="Barry-McIntyre estimator")
   plot(X, add=TRUE, pch=16)
   integral(Z) # should equal 25
}
\keyword{spatial}
\keyword{smooth}
