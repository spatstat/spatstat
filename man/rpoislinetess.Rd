\name{rpoislinetess}
\alias{rpoislinetess}
\title{Poisson Line Tessellation}
\description{
  Generate a tessellation delineated by the lines of
  the Poisson line process
}
\usage{
rpoislinetess(lambda, win = owin())
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson line process.
    A positive number.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
    Currently, the window must be a rectangle.
  }
}
\details{
  This algorithm generates a realisation
  of the uniform Poisson line process, and divides the window
  \code{win} into tiles separated by these lines.

  The argument \code{lambda} must be a positive number.
  It controls the intensity of the process. The expected number of
  lines intersecting a convex region of the plane is equal to
  \code{lambda} times the perimeter length of the region.
  The expected total length of the lines crossing a region of the plane
  is equal to \code{lambda * pi} times the area of the region.
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{rpoisline}}
  to generate the lines only.
}
\examples{
 X <- rpoislinetess(3)
 plot(as.im(X), main="rpoislinetess(3)")
 plot(X, add=TRUE)
}
\keyword{spatial}
\keyword{datagen}
