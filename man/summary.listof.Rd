\name{summary.listof}
\alias{summary.listof}
\title{Summary of a List of Things}
\description{
  Prints a useful summary of each item in a list of things.
}
\usage{
 \method{summary}{listof}(object, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{"listof"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link{summary}}.

  An object of the class \code{"listof"} is effectively a list
  of things which are all of the same class.

  This function extracts a useful summary of each of the items in the list.
}
\seealso{
  \code{\link{summary}},
  \code{\link{plot.listof}}
}
\examples{
  x <- list(A=runif(10), B=runif(10), C=runif(10))
  class(x) <- c("listof", class(x))
  summary(x)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}

