\name{lengths.psp}
\alias{lengths.psp}
\title{Lengths of Line Segments}
\description{
  Computes the length of each line segment
  in a line segment pattern.
}
\usage{
  lengths.psp(x, squared=FALSE)
}
\arguments{
  \item{x}{
    A line segment pattern (object of class \code{"psp"}).
  }
  \item{squared}{
    Logical value indicating whether to return
    the squared lengths (\code{squared=TRUE})
    or the lengths themselves (\code{squared=FALSE}, the default).
  }
}
\value{
  Numeric vector.
}
\details{
  The length of each line segment is computed
  and the lengths are returned as a numeric vector.

  Using squared lengths may be more efficient for some purposes,
  for example, to find the length of the shortest segment,
  \code{sqrt(min(lengths.psp(x, squared=TRUE)))}
  is faster than \code{min(lengths.psp(x))}.
}
\seealso{
  \code{\link{marks.psp}},
  \code{\link{summary.psp}},
  \code{\link{midpoints.psp}},
  \code{\link{angles.psp}},
  \code{\link{endpoints.psp}},
  \code{\link{extrapolate.psp}}.
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  b <- lengths.psp(a)   
}
\author{
  \adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}

