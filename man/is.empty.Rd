\name{is.empty}  
\alias{is.empty}
\alias{is.empty.owin}
\alias{is.empty.ppp}
\alias{is.empty.psp}
\alias{is.empty.default}
\title{Test Whether An Object Is Empty}
\description{
  Checks whether the argument is an empty window,
  an empty point pattern, etc.
}
\usage{
is.empty(x)
\method{is.empty}{owin}(x)
\method{is.empty}{ppp}(x)
\method{is.empty}{psp}(x)
\method{is.empty}{default}(x)
}
\arguments{
  \item{x}{
    A window (object of class \code{"owin"}),
    a point pattern (object of class \code{"ppp"}), or
    a line segment pattern (object of class \code{"psp"}).
  }
}
\details{
  This function tests whether the object \code{x}
  represents an empty spatial object, such as an empty window,
  a point pattern with zero points, or a line segment pattern
  with zero line segments. 

  An empty window can be obtained as the output of
  \code{\link{intersect.owin}}, \code{\link{erosion}},
  \code{\link{opening}}, \code{\link{complement.owin}}
  and some other operations.

  An empty point pattern or line segment pattern can be
  obtained as the result of simulation.
}
\value{
  Logical value.
}

\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
