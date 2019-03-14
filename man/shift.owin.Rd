\name{shift.owin}
\alias{shift.owin}
\title{Apply Vector Translation To Window}
\description{
  Applies a vector shift to a window
}
\usage{
 \method{shift}{owin}(X, vec=c(0,0), \dots, origin=NULL)
}
\arguments{
  \item{X}{Window (object of class \code{"owin"}).}
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
  Another window (of class \code{"owin"}) representing the
  result of applying the vector shift.
}
\details{
  The window is translated by the vector \code{vec}.
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
  \code{\link{shift.ppp}},
  \code{\link{periodify}},
  \code{\link{rotate}},
  \code{\link{affine}},
  \code{\link{centroid.owin}}
}
\examples{
  W <- owin(c(0,1),c(0,1))
  X <- shift(W, c(2,3))
  \dontrun{
  plot(W)
  # no discernible difference except coordinates are different
  }
  shift(W, origin="top")
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
