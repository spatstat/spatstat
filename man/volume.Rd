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
  (two-dimensional windows), which is identical to
  \code{\link{area.owin}},
  and a method for the class \code{"linnet"} of linear networks,
  which returns the length of the network.
}
\seealso{
  \code{\link{area.owin}},
  \code{\link{volume.box3}},
  \code{\link{volume.boxx}},
  \code{\link{volume.linnet}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{math}
