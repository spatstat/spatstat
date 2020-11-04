\name{Extract.ppx}
\alias{[.ppx}
\title{Extract Subset of Multidimensional Point Pattern}
\description{
  Extract a subset of a multidimensional point pattern.
}
\usage{
  \method{[}{ppx}(x, i, drop=FALSE, clip=FALSE, ...)
}
\arguments{
  \item{x}{
    A multidimensional point pattern (object of class \code{"ppx"}).
  }
  \item{i}{
   Subset index. A valid subset index in the usual \R sense,
   indicating which points should be retained;
   or a spatial domain of class \code{"boxx"} or \code{"box3"}.
  }
  \item{drop}{
    Logical value indicating whether to remove unused levels
    of the marks, if the marks are a factor.
  }
  \item{clip}{
    Logical value indicating how to form the domain of the resulting
    point pattern, when \code{i} is a box (object of class \code{"boxx"}). 
    If \code{clip=FALSE} (the default), the result has domain
    equal to \code{i}. If \code{clip=TRUE}, the resulting domain
    is the intersection between the domain of \code{x} and the
    domain \code{i}.
  }
  \item{\dots}{
    Ignored. 
  }
}
\value{
  A multidimensional point pattern (of class \code{"ppx"}).
}
\details{
  This function extracts a designated subset of a multidimensional
  point pattern.

  The function \code{[.ppx} is a method for \code{\link{[}} for the
  class \code{"ppx"}. It extracts a designated subset of a point pattern.
  The argument \code{i} may be either
  \itemize{
    \item a subset index in the usual \R sense:
    either a numeric vector
    of positive indices (identifying the points to be retained),
    a numeric vector of negative indices (identifying the points
    to be deleted) or a logical vector of length equal to the number of
    points in the point pattern \code{x}. In the latter case, 
    the points \code{(x$x[i], x$y[i])} for which 
    \code{subset[i]=TRUE} will be retained, and the others
    will be deleted.
  \item
    a spatial domain of class \code{"boxx"} or \code{"box3"}.
    Points falling inside this region will be retained.
  }
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
   Y <- ppx(coords(cells), domain = boxx(c(0,1),c(0,1)))
   dom <- shift(domain(Y), vec = c(.5,.5))
   Y[dom]
   Y[dom, clip=TRUE]
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
