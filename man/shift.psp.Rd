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
  \item{origin}{Character string determining a location
    that will be shifted to the origin. Options are
    \code{"centroid"}, \code{"midpoint"} and \code{"bottomleft"}.
    Partially matched.
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
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
