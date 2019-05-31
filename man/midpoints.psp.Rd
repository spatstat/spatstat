\name{midpoints.psp}
\alias{midpoints.psp}
\title{Midpoints of Line Segment Pattern}
\description{
  Computes the midpoints of each line segment
  in a line segment pattern.
}
\usage{
  midpoints.psp(x)
}
\arguments{
  \item{x}{
    A line segment pattern (object of class \code{"psp"}).
  }
}
\value{
  Point pattern (object of class \code{"ppp"}).
}
\details{
  The midpoint of each line segment is computed.
}
\seealso{
  \code{\link{marks.psp}},
  \code{\link{summary.psp}},
  \code{\link{lengths.psp}}
  \code{\link{angles.psp}},
  \code{\link{endpoints.psp}},
  \code{\link{extrapolate.psp}}.
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  b <- midpoints.psp(a)   
}
\author{
  \adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}

