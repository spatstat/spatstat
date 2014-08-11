\name{selfcrossing.psp}
\alias{selfcrossing.psp}
\title{Crossing Points in a Line Segment Pattern}
\description{
  Finds any crossing points between the line segments in a 
  line segment pattern.
}
\usage{
  selfcrossing.psp(A)
}
\arguments{
  \item{A}{
    Line segment pattern (object of class \code{"psp"}).
  }
}
\value{
  Point pattern (object of class \code{"ppp"}).
}
\details{
  This function finds any crossing points between
  different line segments in the line segment pattern \code{A}.

  A crossing point occurs whenever one of the line segments in \code{A}
  intersects another line segment in \code{A}, at a nonzero
  angle of intersection.
}
\seealso{
  \code{\link{crossing.psp}},
  \code{\link{psp.object}},
  \code{\link{ppp.object}}.
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  plot(a, col="green", main="selfcrossing.psp")
  P <- selfcrossing.psp(a)
  plot(P, add=TRUE, col="red")
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}

