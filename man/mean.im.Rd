\name{mean.im}
\alias{max.im}
\alias{min.im}
\alias{mean.im}
\alias{median.im}
\alias{range.im}
\alias{sum.im}
\title{Maximum, Minimum, Mean, Median, Range or Sum
  of Pixel Values in an Image}
\description{
  Calculates the mean, median, range, sum, maximum or minimum
  of the pixel values in a pixel image.
}
\usage{
  \method{max}{im}(x, \dots)
  \method{min}{im}(x, \dots)
  \method{mean}{im}(x, \dots)
  \method{median}{im}(x, \dots)
  \method{range}{im}(x, \dots)
  \method{sum}{im}(x, \dots)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"im"}).}
  \item{\dots}{Arguments passed to \code{\link{mean.default}}.}
}
\details{
  These functions calculate the mean, median, range, sum, maximum
  or minimum of the pixel values in the image
  \code{x}.

  An object of class \code{"im"}
  describes a pixel image. See \code{\link{im.object}})
  for details of this class.

  The function \code{mean.im} is a method for the generic
  function \code{\link{mean}} for the class \code{"im"}.
  Similarly \code{median.im} is a method for the generic
  \code{\link{median}} and \code{range.im} is a method for \code{\link{range}}.
  
  If the image \code{x} is logical-valued, the mean value of \code{x} is
  the fraction of pixels that have the value \code{TRUE}. The median is
  not defined. 

  If the image \code{x} is factor-valued, then the mean of \code{x}
  is the mean of the integer codes of the pixel values. The median and range
  are not defined.

  Any arguments in \code{...} are passed to the default method,
  for example \code{\link{mean.default}}.
  In particular, using the argument \code{trim} will compute the
  trimmed mean, as explained in the help for \code{\link{mean.default}}.

  Other information about an image can be obtained using
  \code{\link{summary.im}} or \code{\link{quantile.im}}.
}
\value{
  A single number.
}
\seealso{
  \code{\link{mean}},
  \code{\link{median}},
  \code{\link{range}},
  \code{\link{sum}},
  \code{\link{mean.default}},
  \code{\link{median.default}},
  \code{\link{range.default}},
  \code{\link{quantile.im}},
  \code{\link{im.object}},
  \code{\link{summary.im}}.
}
\examples{
  X <- as.im(function(x,y) {x^2}, unit.square())
  mean(X)
  median(X)
  range(X)

  mean(X, trim=0.05)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
