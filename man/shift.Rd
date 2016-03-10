\name{shift}
\alias{shift}
\title{Apply Vector Translation}
\description{
  Applies a vector shift of the plane 
  to a geometrical object,
  such as a point pattern or a window. 
}
\usage{
  shift(X, \dots)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).}
  \item{\dots}{Arguments determining the shift vector.}
}
\value{
  Another object of the same type, representing the
  result of applying the shift.
}
\details{
  This is generic. Methods are provided for
  point patterns (\code{\link{shift.ppp}})
  and windows (\code{\link{shift.owin}}).

  The object is translated by the vector \code{vec}.
}
\seealso{
  \code{\link{shift.ppp}},
  \code{\link{shift.owin}},
  \code{\link{rotate}},
  \code{\link{affine}},
  \code{\link{periodify}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
