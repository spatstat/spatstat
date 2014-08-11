\name{diameter}
\alias{diameter}
\title{Diameter of an Object}
\description{
  Computes the diameter of an object such as a two-dimensional window
  or three-dimensional box. 
}
\usage{
 diameter(x)
}
\arguments{
  \item{x}{
    A window or other object whose diameter will be computed.
  }
}
\value{
  The numerical value of the diameter of the object.
}
\details{
  This function computes the diameter of an object
  such as a two-dimensional window or a three-dimensional box. 
  The diameter is the maximum distance 
  between any two points in the object.

  The function \code{diameter} is generic, with methods for
  the class \code{"owin"} (two-dimensional windows),
  \code{"box3"} (three-dimensional boxes) and
  \code{"boxx"} (multi-dimensional boxes).
}
\seealso{
  \code{\link{diameter.owin}},
  \code{\link{diameter.box3}},
  \code{\link{diameter.boxx}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
