\name{Extract.ppx}
\alias{[.ppx}
\title{Extract Subset of Multidimensional Point Pattern}
\description{
  Extract a subset of a multidimensional point pattern.
}
\usage{
  \method{[}{ppx}(x, i, drop=FALSE, ...)
}
\arguments{
  \item{x}{
    A multidimensional point pattern (object of class \code{"ppx"}).
  }
  \item{i}{
   Subset index. A valid subset index in the usual \R sense,
   indicating which points should be retained.
  }
  \item{drop}{
    Logical value indicating whether to remove unused levels
    of the marks, if the marks are a factor.
  }
  \item{\dots}{
    Ignored. 
  }
}
\value{
  A multidimensional point pattern (of class \code{"ppx"})
  in the same domain.
}
\details{
  This function extracts a designated subset of a multidimensional
  point pattern.

  The function \code{[.ppx} is a method for \code{\link{[}} for the
  class \code{"ppx"}. It extracts a designated subset of a point pattern.
  The argument \code{i} should be a subset index in the usual \R sense:
  either a numeric vector
  of positive indices (identifying the points to be retained),
  a numeric vector of negative indices (identifying the points
  to be deleted) or a logical vector of length equal to the number of
  points in the point pattern \code{x}. In the latter case, 
  the points \code{(x$x[i], x$y[i])} for which 
  \code{subset[i]=TRUE} will be retained, and the others
  will be deleted.
 
  The argument \code{drop} determines whether to remove
  unused levels of a factor, if the point pattern is multitype
  (i.e. the marks are a factor) or if the marks are a data frame or hyperframe
  in which some of the columns are factors.

  Use the function \code{\link{unmark}} to remove marks from a
  marked point pattern.
}
\seealso{
  \code{\link{ppx}}
}
\examples{
   df <- data.frame(x=runif(4),y=runif(4),z=runif(4))
   X <- ppx(data=df, coord.type=c("s","s","t"))
   X[-2]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
