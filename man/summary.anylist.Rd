\name{summary.anylist}
\alias{summary.anylist}
\title{Summary of a List of Things}
\description{
  Prints a useful summary of each item in a list of things.
}
\usage{
 \method{summary}{anylist}(object, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{"anylist"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link{summary}}.

  An object of the class \code{"anylist"} is effectively a list
  of things which are intended to be treated in a similar way.
  See \code{\link{anylist}}.

  This function extracts a useful summary of each of the items in the list.
}
\seealso{
  \code{\link{anylist}},
  \code{\link{summary}},
  \code{\link{plot.anylist}}
}
\examples{
  x <- anylist(A=runif(10), B=runif(10), C=runif(10))
  summary(x)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{methods}

