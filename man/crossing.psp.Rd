\name{crossing.psp}
\alias{crossing.psp}
\title{Crossing Points of Two Line Segment Patterns}
\description{
  Finds any crossing points between 
  two line segment patterns.
}
\usage{
  crossing.psp(A,B,fatal=TRUE,details=FALSE)
}
\arguments{
  \item{A,B}{
    Line segment patterns (objects of class \code{"psp"}).
  }
  \item{details}{
    Logical value indicating whether to return additional information.
    See below.
  }
  \item{fatal}{
    Logical value indicating what to do
    if the windows of \code{A} and \code{B} do not overlap.
    See Details.
  }
}
\value{
  Point pattern (object of class \code{"ppp"}).
}
\details{
  This function finds any crossing points between
  the line segment patterns \code{A} and \code{B}.

  A crossing point occurs whenever one of the line segments in \code{A}
  intersects one of the line segments in \code{B}, at a nonzero
  angle of intersection.

  The result is a point pattern consisting of all the intersection points.

  If \code{details=TRUE}, additional information is computed,
  specifying where each intersection point came from.
  The resulting point pattern has a data frame of marks, with columns
  named \code{iA, jB, tA, tB}. The marks \code{iA}
  and \code{jB} are the indices of the line segments in \code{A} and
  \code{B}, respectively, which produced each intersection point.
  The marks \code{tA} and \code{tB} are numbers between 0 and 1
  specifying the position of the intersection point along the
  original segments.

  If the windows \code{Window(A)} and \code{Window(B)} do not
  overlap, then an error will be reported if \code{fatal=TRUE},
  while if \code{fatal=FALSE} an error will not occur
  and the result will be \code{NULL}.
}
\seealso{
  \code{\link{selfcrossing.psp}},
  \code{\link{psp.object}},
  \code{\link{ppp.object}}.
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  b <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  plot(a, col="green", main="crossing.psp")
  plot(b, add=TRUE, col="blue")
  P <- crossing.psp(a,b)
  plot(P, add=TRUE, col="red")
  as.data.frame(crossing.psp(a,b,details=TRUE))
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{manip}

