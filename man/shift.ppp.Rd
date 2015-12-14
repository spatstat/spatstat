\name{shift.ppp}
\alias{shift.ppp}
\title{Apply Vector Translation To Point Pattern}
\description{
  Applies a vector shift to a point pattern. 
}
\usage{
 \method{shift}{ppp}(X, vec=c(0,0), \dots, origin=NULL)
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"}).}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
  \item{origin}{Character string determining a location
    that will be shifted to the origin. Options are
    \code{"centroid"}, \code{"midpoint"} and \code{"bottomleft"}.
    Partially matched.
  }
}
\value{
  Another point pattern (of class \code{"ppp"}) representing the
  result of applying the vector shift.
}
\details{
  The point pattern, and its window, are
  translated by the vector \code{vec}.
  
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
  \code{\link{shift.owin}},
  \code{\link{periodify}},
  \code{\link{rotate}},
  \code{\link{affine}}
}
\examples{
  data(cells)
  X <- shift(cells, c(2,3))
  \dontrun{
  plot(X)
  # no discernible difference except coordinates are different
  }
  plot(cells, pch=16)
  plot(shift(cells, c(0.03,0.03)), add=TRUE)

  shift(cells, origin="mid")
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
