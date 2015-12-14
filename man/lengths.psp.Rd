\name{lengths.psp}
\alias{lengths.psp}
\title{Lengths of Line Segments}
\description{
  Computes the length of each line segment
  in a line segment pattern.
}
\usage{
  lengths.psp(x)
}
\arguments{
  \item{x}{
    A line segment pattern (object of class \code{"psp"}).
  }
}
\value{
  Numeric vector.
}
\details{
  The length of each line segment is computed
  and the lengths are returned as a numeric vector.
}
\seealso{
  \code{\link{summary.psp}},
  \code{\link{midpoints.psp}},
  \code{\link{angles.psp}}
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  b <- lengths.psp(a)   
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}

