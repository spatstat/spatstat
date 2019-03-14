\name{shift.psp}
\alias{shift.psp}
\title{Apply Vector Translation To Line Segment Pattern}
\description{
  Applies a vector shift to a line segment pattern. 
}
\usage{
 \method{shift}{psp}(X, vec=c(0,0), \dots, origin=NULL)
}
\arguments{
  \item{X}{Line Segment pattern (object of class \code{"psp"}).}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
  \item{origin}{
    Location that will be shifted to the origin.
    Either a numeric vector of length 2 giving the location,
    or a point pattern containing only one point,
    or a list with two entries named \code{x} and \code{y},
    or one of the character strings
    \code{"centroid"}, \code{"midpoint"},
    \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"},
    \code{"topleft"}, \code{"bottomleft"}, \code{"topright"} or
    \code{"bottomright"} (partially matched).
  }
}
\value{
  Another line segment pattern (of class \code{"psp"}) representing the
  result of applying the vector shift.
}
\details{
  The line segment pattern, and its window, are
  translated by the vector \code{vec}.
  
  This is a method for the generic function \code{\link{shift}}.

  If \code{origin} is given,
  the argument \code{vec} will be ignored; instead the shift will be performed
  so that the specified geometric location is shifted to the
  coordinate origin \eqn{(0,0)}.
  The argument \code{origin} should be either a numeric vector of length
  2 giving the spatial coordinates of a location, or one of the character
  strings \code{"centroid"}, \code{"midpoint"},
  \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"},
  \code{"topleft"}, \code{"bottomleft"}, \code{"topright"} or
  \code{"bottomright"} (partially matched).
  If \code{origin="centroid"} then the centroid of the window will be
  shifted to the origin. If \code{origin="midpoint"} then the centre of
  the bounding rectangle of the window will be shifted to the origin.
  If \code{origin="bottomleft"} then the bottom left corner of the
  bounding rectangle of the window will be shifted to the origin,
  and so on.
}
\seealso{
  \code{\link{shift}},
  \code{\link{shift.owin}},
  \code{\link{shift.ppp}},
  \code{\link{periodify}},
  \code{\link{rotate}},
  \code{\link{affine}}
}
\examples{
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  plot(X, col="red")
  Y <- shift(X, c(0.05,0.05))
  plot(Y, add=TRUE, col="blue")

  shift(Y, origin="mid")
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
