\name{mean.linim}
\alias{mean.linim}
\alias{median.linim}
\alias{quantile.linim}
\title{Mean, Median, Quantiles of Pixel Values on a Linear Network}
\description{
  Calculates the mean, median, or quantiles
  of the pixel values in a pixel image on a linear network.
}
\usage{
  \method{mean}{linim}(x, \dots)

  \method{median}{linim}(x, \dots)

  \method{quantile}{linim}(x, probs=seq(0,1,0.25), \dots)
}
\arguments{
  \item{x}{
    A pixel image on a linear network (object of class
    \code{"linim"}).
  }
  \item{probs}{
    Vector of probabilities for which quantiles should be
    calculated.
  } 
  \item{\dots}{Arguments passed to other methods.}
}
\details{
  These functions calculate the mean, median and quantiles
  of the pixel values in the image
  \code{x} on a linear network.

  An object of class \code{"linim"}
  describes a pixel image on a linear network. See \code{\link{linim}}.

  The functions described here are methods for the 
  generic \code{\link{mean}}, \code{\link[stats]{median}}
  and \code{\link[stats]{quantile}} for the class \code{"linim"}.
}
\value{
  For \code{mean} and \code{median}, a single number.
  For \code{quantile}, a numeric vector of the same length as \code{probs}.
}
\seealso{
  \code{\link{mean}},
  \code{\link[stats]{median}},
  \code{\link[stats]{quantile}},

  \code{\link{mean.im}}.
}
\examples{
  M <- as.mask.psp(as.psp(simplenet))
  Z <- as.im(function(x,y) {x-y}, W=M)
  X <- linim(simplenet, Z)
  X
  mean(X)
  median(X)
  quantile(X)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
