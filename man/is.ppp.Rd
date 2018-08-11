\name{is.ppp}
\alias{is.ppp}
\title{Test Whether An Object Is A Point Pattern}
\description{
  Checks whether its argument is a point pattern
  (object of class \code{"ppp"}).
}
\usage{
is.ppp(x)
}
\arguments{
  \item{x}{Any object.}
}
\details{
  This function tests whether the object \code{x}
  is a point pattern object of class
  \code{"ppp"}. See \code{\link{ppm.object}} for details
  of this class.

  The result is determined to be \code{TRUE} if \code{x}
  inherits from \code{"ppp"}, i.e. if \code{x}
  has \code{"ppp"} amongst its classes.
}
\value{
  \code{TRUE} if \code{x} is a point pattern,
  otherwise \code{FALSE}.
}

\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
