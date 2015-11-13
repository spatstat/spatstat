\name{Extract.solist}
\alias{[.solist}
\alias{[<-.solist}
\title{Extract or Replace Subset of a List of Spatial Objest}
\description{
  Replace a subset of a list of things.
}
\usage{
  \method{[}{solist}(x, i, \dots)

  \method{[}{solist}(x, i) <- value
}
\arguments{
  \item{x}{
    An object of class \code{"solist"} representing a list of
    two-dimensional spatial objects.
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
  Another object of class \code{"solist"}.
}
\details{
  This is a subset replacement method for the class \code{"solist"}.

  The argument \code{x} should be an object of class \code{"solist"}
  representing a list of two-dimensional spatial objects.
  See \code{\link{solist}}.

  The method replaces a designated
  subset of \code{x}, and returns an object of class \code{"solist"}.
}
\seealso{
  \code{\link{solist}},
  \code{\link{plot.solist}},
  \code{\link{summary.solist}}
}
\examples{
   x <- solist(japanesepines, cells, redwood)
   x[1] <- list(finpines)
 }
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{list}
\keyword{manip}
