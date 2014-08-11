\name{append.psp}
\alias{append.psp}
\title{Combine Two Line Segment Patterns}
\description{
  Combine two line segment patterns into a single pattern.
}
\usage{
  append.psp(A, B)
}
\arguments{
  \item{A,B}{
    Line segment patterns (objects of class \code{"psp"}).
  }
}
\value{
  Another line segment pattern (object of class \code{"psp"}).
}
\details{
  This function is used to superimpose two line segment patterns
  \code{A} and \code{B}.
  
  The two patterns must have \bold{identical} windows. If one
  pattern has marks, then the other must also have marks
  of the same type.  It the marks are data frames then the
  number of columns of these data frames, and the names of
  the columns must be identical.

  (To combine two point patterns, see \code{superimpose}).
}
\seealso{
  \code{\link{psp}},
  \code{\link{as.psp}},
  \code{\link{superimpose}},
}
\examples{
  X <- psp(runif(20), runif(20), runif(20), runif(20),  window=owin())
  Y <- psp(runif(5), runif(5), runif(5), runif(5),  window=owin())
  append.psp(X,Y)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
