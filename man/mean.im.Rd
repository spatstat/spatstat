\name{mean.im}
\alias{mean.im}
\alias{median.im}
\title{Mean and Median of Pixel Values in an Image}
\description{
  Calculates the mean or median
  of the pixel values in a pixel image.
}
\usage{
  \method{mean}{im}(x, trim=0, na.rm=TRUE, \dots)

  \method{median}{im}(x, na.rm=TRUE)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"im"}).}
  \item{na.rm}{
    Logical value indicating whether \code{NA} values should be
    stripped before the computation proceeds.
  }
  \item{trim}{
    The fraction (0 to 0.5) of pixel values to be trimmed from
    each end of their range, before the mean is computed.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  These functions calculate the mean and median
  of the pixel values in the image \code{x}.

  An object of class \code{"im"}
  describes a pixel image. See \code{\link{im.object}})
  for details of this class.

  The function \code{mean.im} is a method for the generic
  function \code{\link[base]{mean}} for the class \code{"im"}.
  Similarly \code{median.im} is a method for the generic
  \code{\link[stats]{median}}.
  
  If the image \code{x} is logical-valued, the mean value of \code{x} is
  the fraction of pixels that have the value \code{TRUE}. The median is
  not defined. 

  If the image \code{x} is factor-valued, then the mean of \code{x}
  is the mean of the integer codes of the pixel values. The median is
  are not defined.

  Other mathematical operations on images are supported by
  \code{\link{Math.im}}, \code{\link{Summary.im}}
  and \code{\link{Complex.im}}.

  Other information about an image can be obtained using
  \code{\link{summary.im}} or \code{\link{quantile.im}}.
}
\value{
  A single number.
}
\seealso{
  \code{\link{Math.im}} for other operations.

  Generics and default methods:
  \code{\link[base]{mean}},
  \code{\link[stats]{median}}.
  
  \code{\link{quantile.im}},
  \code{\link{anyNA.im}},
  \code{\link{im.object}},
  \code{\link{summary.im}}.
}
\examples{
  X <- as.im(function(x,y) {x^2}, unit.square())
  mean(X)
  median(X)
  mean(X, trim=0.05)
}
\author{
  \spatstatAuthors and Kassel Hingee.
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
