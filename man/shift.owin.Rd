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
  \item{origin}{Character string determining a location
    that will be shifted to the origin. Options are
    \code{"centroid"}, \code{"midpoint"} and \code{"bottomleft"}.
    Partially matched. 
  }
}
\value{
  Another window (of class \code{"owin"}) representing the
  result of applying the vector shift.
}
\details{
  The window is translated by the vector \code{vec}.
  This is a method for the generic function \code{\link{shift}}.

  If \code{origin} is given, then it should be one of the character
  strings \code{"centroid"}, \code{"midpoint"} or \code{"bottomleft"}.
  The argument \code{vec} will be ignored; instead the shift will be performed
  so that the specified geometric location is shifted to the origin.
  If \code{origin="centroid"} then the centroid of the window will be
  shifted to the origin. If \code{origin="midpoint"} then the centre of
  the bounding rectangle of the window will be shifted to the origin.
  If \code{origin="bottomleft"} then the bottom left corner of the
  bounding rectangle of the window will be shifted to the origin.
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
  shift(W, origin="mid")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
