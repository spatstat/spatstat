\name{vertices}
\alias{vertices}
\alias{vertices.owin}
\title{Vertices of a Window}
\description{
  Finds the vertices of a window, or similar object.
}
\usage{
 vertices(w)

 \method{vertices}{owin}(w)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}) or similar object.}
}
\value{
  A list with components \code{x} and \code{y} giving the coordinates
  of the vertices.
}
\details{
  This function computes the vertices (`corners') of a spatial window
  or other object.

  For \code{vertices.owin},
  the argument \code{w} should be a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details).

  If \code{w} is a rectangle, the coordinates of the
  four corner points are returned.

  If \code{w} is a polygonal window (consisting of one or more
  polygons), the coordinates of the vertices of all polygons are
  returned.

  If \code{w} is a binary mask, then a `boundary pixel' is defined to be
  a pixel inside the window
  which has at least one neighbour outside the window. The coordinates
  of the centres of all boundary pixels are returned.
}
\seealso{
  \code{\link{owin.object}}.
}
\examples{
  data(letterR)
  vert <- vertices(letterR)

  plot(letterR, main="Polygonal vertices")
  points(vert)
  plot(letterR, main="Boundary pixels")
  points(vertices(as.mask(letterR)))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
