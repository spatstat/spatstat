\name{Extract.lpp}
\alias{[.lpp}
\title{Extract Subset of Point Pattern on Linear Network}
\description{
  Extract a subset of a point pattern on a linear network.
}
\usage{
  \method{[}{lpp}(x, i, ...)
}
\arguments{
  \item{x}{
    A point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{i}{
    Subset index. A valid subset index in the usual \R sense,
   indicating which points should be retained.
  }
  \item{\dots}{
    Ignored. 
  }
}
\value{
  A point pattern on a linear network (of class \code{"lpp"}).
}
\details{
  This function extracts a designated subset of a point pattern
  on a linear network.

  The function \code{[.lpp} is a method for \code{\link{[}} for the
  class \code{"lpp"}. It extracts a designated subset of a point pattern.
  The argument \code{i} should be a subset index in the usual \R sense:
  either a numeric vector
  of positive indices (identifying the points to be retained),
  a numeric vector of negative indices (identifying the points
  to be deleted) or a logical vector of length equal to the number of
  points in the point pattern \code{x}. In the latter case, 
  the points \code{(x$x[i], x$y[i])} for which 
  \code{subset[i]=TRUE} will be retained, and the others
  will be deleted.
 
  Use the function \code{\link{unmark}} to remove marks from a
  marked point pattern.
}
\seealso{
  \code{\link{lpp}}
}
\examples{
  # Chicago crimes data - remove cases of assault
  chicago[marks(chicago) != "assault"]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
