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
  \code{\link{summary.psp}},
  \code{\link{lengths.psp}},
  \code{\link{angles.psp}}
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  b <- midpoints.psp(a)   
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}

