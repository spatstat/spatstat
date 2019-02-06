\name{densityAdaptiveKernel}
\alias{densityAdaptiveKernel}
\alias{densityAdaptiveKernel.ppp}
\title{Adaptive Kernel Estimate of Intensity of Point Pattern}
\description{
  Computes an adaptive estimate of the intensity function of a point
  pattern using a variable-bandwidth smoothing kernel.
}
\usage{
densityAdaptiveKernel(X, \dots)

\method{densityAdaptiveKernel}{ppp}(X, \dots)
}
\arguments{
  \item{X}{Point pattern dataset (object of class \code{"ppp"}).}
  \item{\dots}{Arguments passed to \code{\link{as.im}} determining the
    pixel resolution of the result.
  }
}
\details{
  This function is an alternative to \code{\link{density.ppp}}. It
  computes an estimate of the intensity function of a point pattern
  dataset. The result is a pixel image giving the estimated intensity.

  \bold{This is currently a stub which just calls} \code{\link{density.ppp}}.
}
\value{
  A pixel image (object of class \code{"im"}) whose values are
  estimates of the intensity of \code{X}.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{adaptive.density}},
  \code{\link{densityVoronoi}},
  \code{\link{im.object}}.
}
\examples{
  plot(densityAdaptiveKernel(nztrees), main="Adaptive kernel estimate")
}
\author{
  \spatstatAuthors and Tilman Davies
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
