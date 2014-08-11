\name{simplify.owin}
\Rdversion{1.1}
\alias{simplify.owin}
\title{
  Approximate a Polygon by a Simpler Polygon
}
\description{
  Given a polygonal window, this function finds a simpler polygon
  that approximates it.
}
\usage{
simplify.owin(W, dmin)
}
\arguments{
  \item{W}{
    The polygon which is to be simplied.
    An object of class \code{"owin"}.
  }
  \item{dmin}{
    Numeric value. The smallest permissible length of an edge.
  }
}
\details{
  This function simplifies a polygon \code{W}
  by recursively deleting the shortest edge of \code{W}
  until all remaining edges are longer than the specified
  minimum length \code{dmin}, or until there are only three edges left.

  The argument \code{W} must be a window (object of class
  \code{"owin"}). It should be of type \code{"polygonal"}.
  If \code{W} is a rectangle, it is returned without alteration.

  The simplification algorithm is not yet implemented for
  binary masks. If \code{W} is a mask, an error is generated.
}
\value{
  Another window (object of class \code{"owin"})
  of type \code{"polygonal"}.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{owin}}
}
\examples{
  data(letterR)
  plot(letterR, col="red")
  plot(simplify.owin(letterR, 0.3), col="blue", add=TRUE)
  data(chorley)
  W <- chorley$window
  plot(W)
  WS <- simplify.owin(W, 2)
  plot(WS, add=TRUE, border="green")
  points(vertices(WS))
}
\keyword{spatial}
\keyword{math}
