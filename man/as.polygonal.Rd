\name{as.polygonal}
\Rdversion{1.1}
\alias{as.polygonal}
\title{
Convert a Window to a Polygonal Window
}
\description{
  Given a window \code{W} of any geometric type (rectangular,
  polygonal or binary mask),  this function returns a polygonal window
  that represents the same spatial domain.
}
\usage{
as.polygonal(W)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}).
  }
}
\details{
  Given a window \code{W} of any geometric type (rectangular,
  polygonal or binary mask),  this function returns a polygonal window
  that represents the same spatial domain.

  If \code{W} is already polygonal, it is returned without change.

  If \code{W} is a rectangle, it is converted to a polygon with
  4 vertices.

  If \code{W} is a binary mask, then each pixel in the
  mask is replaced by a small square or rectangle, and the union of these
  squares or rectangles is computed. The result is a polygonal window
  that has only horizontal and vertical edges. 
  (Use \code{\link{simplify.owin}} to remove the staircase
  appearance, if desired).
}
\value{
  A polygonal window (object of class \code{"owin"}
  and of type \code{"polygonal"}). 
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{as.mask}},
  \code{\link{simplify.owin}}
}
\examples{
   data(letterR)
   m <- as.mask(letterR, dimyx=32)
   p <- as.polygonal(m)
   if(interactive()) {
      plot(m)
      plot(p, add=TRUE, lwd=2)
   }
}
\keyword{spatial}
\keyword{manip}
