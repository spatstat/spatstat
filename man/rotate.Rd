\name{rotate}
\alias{rotate}
\title{Rotate}
\description{
  Applies a rotation to any two-dimensional object,
  such as a point pattern or a window. 
}
\usage{
  rotate(X, \dots)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).}
  \item{\dots}{Data specifying the rotation.}
}
\value{
  Another object of the same type, representing the
  result of rotating \code{X} through the specified angle.
}
\details{
  This is generic. Methods are provided for
  point patterns (\code{\link{rotate.ppp}})
  and windows (\code{\link{rotate.owin}}).
}
\seealso{
  \code{\link{rotate.ppp}},
  \code{\link{rotate.owin}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
