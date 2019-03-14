\name{shift.im}
\alias{shift.im}
\title{Apply Vector Translation To Pixel Image}
\description{
  Applies a vector shift to a pixel image
}
\usage{
 \method{shift}{im}(X, vec=c(0,0), \dots, origin=NULL)
}
\arguments{
  \item{X}{Pixel image (object of class \code{"im"}).}
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
  Another pixel image (of class \code{"im"}) representing the
  result of applying the vector shift.
}
\details{
  The spatial location of each pixel in the image
  is translated by the vector \code{vec}.
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
  \code{\link{shift}}
}
\examples{
 # make up an image
 X <- setcov(unit.square())
 plot(X)

 Y <- shift(X, c(10,10))
 plot(Y)
 # no discernible difference except coordinates are different

 shift(X, origin="c")
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
