\name{crossing.psp}
\alias{crossing.psp}
\title{Crossing Points of Two Line Segment Patterns}
\description{
  Finds any crossing points between 
  two line segment patterns.
}
\usage{
  crossing.psp(A,B,fatal=TRUE)
}
\arguments{
  \item{A,B}{
    Line segment patterns (objects of class \code{"psp"}).
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
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}

