\name{lengths_psp}
\alias{lengths.psp}
\alias{lengths_psp}
\title{Lengths of Line Segments}
\description{
  Computes the length of each line segment
  in a line segment pattern.
}
\usage{
  lengths.psp(x, squared=FALSE)
  lengths_psp(x, squared=FALSE)
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

  The functions \code{lengths.psp} and \code{lengths_psp} are
  identical. We recommend using the newer name \code{lengths_psp}.
}
\section{Change of name}{
  The name of this function has changed from \code{lengths.psp}
  to \code{lengths_psp}, because the old name \code{lengths.psp}
  could be misinterpreted as a method for \code{\link[base]{lengths}}.
  
  The older function name \code{lengths.psp} is retained temporarily,
  for consistency with older code and documentation.
  
  In future versions of \pkg{spatstat}, the 
  function name \code{lengths.psp} will be removed.
  The newer function name \code{lengths_psp} should be used.
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
  b <- lengths_psp(a)   
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}

