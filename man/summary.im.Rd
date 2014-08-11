\name{summary.im}
\alias{summary.im}
\alias{print.summary.im}
\title{Summarizing a Pixel Image}
\description{
  \code{summary} method for class \code{"im"}.
}
\usage{
  \method{summary}{im}(object, \dots)
  \method{print}{summary.im}(x, \dots)
}
\arguments{
  \item{object}{A pixel image.}
  \item{\dots}{Ignored.}
  \item{x}{Object of class \code{"summary.im"} as returned by
    \code{summary.im}.
  }
}
\details{
  This is a method for the generic \code{\link{summary}}
  for the class \code{"im"}. An object of class \code{"im"}
  describes a pixel image. See \code{\link{im.object}})
  for details of this class. 
  
  \code{summary.im} extracts information about the pixel image,
  and \code{print.summary.im} prints this information in a
  comprehensible format.

  In normal usage, \code{print.summary.im} is invoked implicitly
  when the user calls \code{summary.im} without assigning its value
  to anything. See the examples.

  The information extracted by \code{summary.im} includes
  \describe{
    \item{range}{The range of the image values.}
    \item{mean}{The mean of the image values.}
    \item{integral}{The ``integral'' of the image values,
      calculated as the sum of the image values
      multiplied by the area of one pixel.}
    \item{dim}{The dimensions of the pixel array:
      \code{dim[1]} is the number of rows in the array,
      corresponding to the \bold{y} coordinate.}
  }
}
\value{
  \code{summary.im} returns an object of class \code{"summary.im"},
  while \code{print.summary.im} returns \code{NULL}.
}
\examples{
  # make an image
  X <- as.im(function(x,y) {x^2}, unit.square())
  # summarize it
  summary(X)
  # save the summary
  s <- summary(X)
  # print it
  print(X)
  s
  # extract stuff
  X$dim
  X$range
  X$integral
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
