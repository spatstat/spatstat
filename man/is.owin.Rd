\name{is.owin}
\alias{is.owin}
\title{Test Whether An Object Is A Window}
\description{
  Checks whether its argument is a window
  (object of class \code{"owin"}).
}
\usage{
is.owin(x)
}
\arguments{
  \item{x}{Any object.}
}
\details{
  This function tests whether the object \code{x}
  is a window object of class
  \code{"owin"}. See \code{\link{owin.object}} for details
  of this class.

  The result is determined to be \code{TRUE} if \code{x}
  inherits from \code{"owin"}, i.e. if \code{x}
  has \code{"owin"} amongst its classes.
}
\value{
  \code{TRUE} if \code{x} is a point pattern,
  otherwise \code{FALSE}.
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
