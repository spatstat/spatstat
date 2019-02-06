\name{adaptive.density}
\alias{adaptive.density}
\title{Adaptive Estimate of Intensity of Point Pattern}
\description{
  Computes an adaptive estimate of the intensity function of a point
  pattern.
}
\usage{
adaptive.density(X, \dots, method=c("voronoi","kernel"))
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"} or
    \code{"lpp"}).}
  \item{method}{Character string specifying the estimation method}
  \item{\dots}{
    Additional arguments passed to
    \code{\link{densityVoronoi}} or \code{\link{densityAdaptiveKernel}}.
  }
}
\details{
  This function is an alternative to \code{\link{density.ppp}}. It
  computes an estimate of the intensity function of a point pattern
  dataset. The result is a pixel image giving the estimated intensity.

  If \code{method="voronoi"} the data are passed to the function
  \code{\link{densityVoronoi}} which estimates the intensity using
  the Voronoi-Dirichlet tessellation.

  If \code{method="kernel"} the data are passed to the function
  \code{\link{densityAdaptiveKernel}} which estimates the intensity
  using a variable-bandwidth kernel estimator.
}
\value{
  A pixel image (object of class \code{"im"}) whose values are
  estimates of the intensity of \code{X}.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{densityVoronoi}},
  \code{\link{densityAdaptiveKernel}},
  \code{\link{im.object}}.
}
\examples{
  plot(adaptive.density(nztrees, 1), main="Voronoi estimate")
}
\author{
  \spatstatAuthors and Mehdi Moradi.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}

    