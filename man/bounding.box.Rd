\name{bounding.box}
\alias{bounding.box}
\title{Bounding Box of a Window or Point Pattern}
\description{
  Find the smallest rectangle containing a given window(s)
  or point pattern(s).
}
\usage{
 bounding.box(\dots)
}
\arguments{
  \item{\dots}{One or more windows (objects of class \code{"owin"}),
    pixel images (objects of class \code{"im"}) or
    point patterns (objects of class \code{"ppp"}).
  }
}
\value{
  A window (object of class \code{"owin"})
  of type \code{"rectangle"} representing a rectangle.
}
\details{
  This function finds the smallest rectangle (with sides parallel to
  the coordinate axes) that contains all the given objects.

  For a window (object of class \code{"owin"}), the bounding box
  is the smallest rectangle that contains all the vertices of the
  window (this is generally smaller than the enclosing frame,
  which is returned by \code{\link{as.rectangle}}).

  For a point pattern (object of class \code{"ppp"}), the bounding box
  is the smallest rectangle that contains all the points of the pattern,
  and is computed by \code{\link{bounding.box.xy}}.

  For a pixel image (object of class \code{"im"}), the image will
  be converted to a window using \code{\link{as.owin}},
  and the bounding box of this window is obtained.
  
  If the argument is a list of several objects, then
  this function finds the smallest rectangle that contains
  all the bounding boxes of the objects.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{as.rectangle}}
  \code{\link{bounding.box.xy}}
}
\examples{
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  r <- bounding.box(w)
  # returns rectangle [1,3] x [2,7]

  w2 <- unit.square()
  r <- bounding.box(w, w2)
  # returns rectangle [0,3] x [0,7]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{utilities}
