\name{compatible}
\alias{compatible}
\title{Test Whether Objects Are Compatible}
\description{
  Tests whether two or more objects of the same class 
  are compatible.
}
\usage{
  compatible(A, B, \dots)
}
\arguments{
  \item{A,B,\dots}{Two or more objects of the same class}
}
\details{
  This generic function is used to check whether the
  objects \code{A} and \code{B} (and any additional
  objects \code{\dots}) are compatible.

  What is meant by \sQuote{compatible} depends on the class of object.

  There are methods for the classes \code{"fv"}, \code{"fasp"},
  \code{"im"} and \code{"units"}.
}
\value{
  Logical value: \code{TRUE} if the objects are compatible, and \code{FALSE}
  if they are not.
}
\seealso{
  \code{\link{compatible.fv}},
  \code{\link{compatible.fasp}},
  \code{\link{compatible.im}},
  \code{\link{compatible.units}} 
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
