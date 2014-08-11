\name{Extract.splitppp}
\alias{[.splitppp}
\alias{[<-.splitppp}
\title{Extract or Replace Sub-Patterns}
\description{
  Extract or replace some of the sub-patterns in a split point pattern.
}
\usage{
  \method{[}{splitppp}(x, ...)
  \method{[}{splitppp}(x, ...) <- value
}
\arguments{
  \item{x}{
    An object of class \code{"splitppp"}, representing a point pattern
    separated into a list of sub-patterns.
  }
  \item{\dots}{
    Subset index. Any valid subset index in the usual \R sense.
  }
  \item{value}{
    Replacement value for the subset. A list of point patterns.
  }
}
\value{
  Another object of class \code{"splitppp"}.
}
\details{
  These are subset methods for the class \code{"splitppp"}.

  The argument \code{x} should be an object of class \code{"splitppp"},
  representing a point pattern that has been separated into a
  list of sub-patterns. It is created by \code{\link{split.ppp}}.
  
  The methods extract or replace a designated
  subset of the list \code{x}, and return an object of class \code{"splitppp"}.
}
\seealso{
  \code{\link{split.ppp}},
  \code{\link{plot.splitppp}},
  \code{\link{summary.splitppp}}
}
\examples{
  data(amacrine)  # multitype point pattern
  y <- split(amacrine)
  y[1]
  y["off"]
  y[1] <- list(runifpoint(42, amacrine$window))
 }
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
