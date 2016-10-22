\name{Smooth.ssf}
\alias{Smooth.ssf}
\title{
  Smooth a Spatially Sampled Function
}
\description{
  Applies kernel smoothing to a spatially sampled function.
}
\usage{
 \method{Smooth}{ssf}(X, \dots)
}
\arguments{
  \item{X}{
    Object of class \code{"ssf"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{Smooth.ppp}}
    to control the smoothing.
  }
}
\details{
  An object of class \code{"ssf"} represents a real-valued or
  vector-valued function that has been evaluated or sampled at an
  irregular set of points.

  The function values will be smoothed using a Gaussian
  kernel.
}
\value{
  A pixel image or a list of pixel images.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{ssf}},
  \code{\link{Smooth.ppp}}
}
\examples{
   f <- ssf(redwood, nndist(redwood))
   Smooth(f, sigma=0.1)
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
