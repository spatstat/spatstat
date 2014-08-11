\name{integral.im}
\alias{integral.im}
\title{
  Integral of a Pixel Image
}
\description{
  Computes the integral of a pixel image.
}
\usage{
integral.im(x, ...)
}
\arguments{
  \item{x}{
    A pixel image (object of class \code{"im"}) with pixel values
    that can be treated as numeric or complex values.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This function treats the pixel image \code{x} as a function of
  the spatial coordinates, and computes its integral.
  The integral is calculated
  by summing the pixel values and multiplying by the area of one pixel.

  The pixel values of \code{x} may be numeric, integer, logical or
  complex. They cannot be factor or character values.

  The logical values \code{TRUE} and \code{FALSE} are converted to
  \code{1} and \code{0} respectively, so that the integral of a logical
  image is the total area of the \code{TRUE} pixels, in the same units
  as \code{unitname(x)}.

  For more complicated integration tasks such as computing the integral
  of an image over a specified subset, use \code{\link{eval.im}} to
  construct an integrand or \code{\link{[.im}} to extract a subset of
  the image.
}
\value{
  A single numeric or complex value.
}
\seealso{
  \code{\link{eval.im}},
  \code{\link{[.im}}
}
\examples{
   # approximate integral of f(x,y) dx dy
   f <- function(x,y){3*x^2 + 2*y}
   Z <- as.im(f, square(1))
   integral.im(Z)
   # correct answer is 2

   data(cells)
   D <- density(cells)
   integral.im(D)
   # should be approximately equal to number of points = 42

   # integrate over the subset [0.1,0.9] x [0.2,0.8]
   W <- owin(c(0.1,0.9), c(0.2,0.8))
   DW <- D[W, drop=FALSE]
   integral.im(DW)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
