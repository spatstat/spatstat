\name{Extract.anylist}
\alias{[.anylist}
\alias{[<-.anylist}
\title{Extract or Replace Subset of a List of Things}
\description{
  Extract or replace a subset of a list of things.
}
\usage{
  \method{[}{anylist}(x, i, \dots)

  \method{[}{anylist}(x, i) <- value
}
\arguments{
  \item{x}{
    An object of class \code{"anylist"} representing a list of things.
  }
  \item{i}{
    Subset index. Any valid subset index in the usual \R sense.
  }
  \item{value}{
    Replacement value for the subset. 
  }
  \item{\dots}{Ignored.}
}
\value{
  Another object of class \code{"anylist"}.
}
\details{
  These are the methods for extracting and replacing subsets 
  for the class \code{"anylist"}.

  The argument \code{x} should be an object of class \code{"anylist"}
  representing a list of things. See \code{\link{anylist}}.

  The method replaces a designated
  subset of \code{x}, and returns an object of class \code{"anylist"}.
}
\seealso{
  \code{\link{anylist}},
  \code{\link{plot.anylist}},
  \code{\link{summary.anylist}}
}
\examples{
   x <- anylist(A=runif(10), B=runif(10), C=runif(10))
   x[1] <- list(A=rnorm(10))
 }
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{list}
\keyword{manip}
