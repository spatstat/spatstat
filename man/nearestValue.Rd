\name{nearestValue}
\alias{nearestValue}
\title{
  Image of Nearest Defined Pixel Value
}
\description{
  Given a pixel image defined on a subset of a rectangle,
  this function assigns a value to every pixel in the rectangle,
  by looking up the value of the nearest pixel that has a value.
}
\usage{
nearestValue(X)
}
\arguments{
  \item{X}{
    A pixel image (object of class \code{"im"}).
  }
}
\details{
  A pixel image in \pkg{spatstat} is always stored on
  a rectangular grid of pixels, but its value may be \code{NA}
  on some pixels, indicating that the image is not defined at those
  pixels.

  This function assigns a value to every pixel in the rectangular grid.
  For each pixel \code{a} in the grid, if the value of \code{X} is not
  defined at \code{a}, the function finds the nearest other pixel \code{b}
  at which the value of \code{X} is defined, and takes the pixel value at
  \code{b} as the new pixel value at \code{a}.
}
\value{
  Another image of the same kind as \code{X}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{blur}}, \code{\link{Smooth.ppp}}
}
\examples{
  X <- as.im(function(x,y) { x + y }, letterR)
  Y <- nearestValue(X)
  plot(solist("X"=X,"nearestValue(X)"=Y), main="", panel.end=letterR)
}
\keyword{spatial}
\keyword{manip}
