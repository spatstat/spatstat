\name{volume}
\alias{volume}
\title{Volume of an Object}
\description{
  Computes the volume of a spatial object such as a three-dimensional box. 
}
\usage{
 volume(x)
}
\arguments{
  \item{x}{
    An object whose volume will be computed.
  }
}
\value{
  The numerical value of the volume of the object.
}
\details{
  This function computes the volume of an object
  such as a three-dimensional box. 

  The function \code{volume} is generic, with methods for
  the classes \code{"box3"} (three-dimensional boxes) and
  \code{"boxx"} (multi-dimensional boxes).
  There is also a method for the class \code{"owin"}
  (two-dimensional windows), which is identical to \code{\link{area.owin}}.

}
\seealso{
  \code{\link{area.owin}},
  \code{\link{volume.box3}},
  \code{\link{volume.boxx}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
