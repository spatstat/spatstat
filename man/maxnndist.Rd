\name{maxnndist}
\alias{maxnndist}
\alias{minnndist}
\title{
  Compute Minimum or Maximum Nearest-Neighbour Distance
}
\description{
  A faster way to compute the minimum or maximum
  nearest-neighbour distance in a point pattern.
}
\usage{
minnndist(X, positive=FALSE, by=NULL)
maxnndist(X, positive=FALSE, by=NULL)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{positive}{
    Logical. If \code{FALSE} (the default), compute the usual
    nearest-neighbour distance.
    If \code{TRUE}, ignore coincident points, so that the
    nearest neighbour distance for each point is greater than zero.
  }
  \item{by}{
    Optional. A factor, which separates \code{X} into groups.
    The algorithm will compute the distance to
    the nearest point in each group. 
  }
}
\details{
  These functions find the minimum and maximum values
  of nearest-neighbour distances in the point pattern \code{X}.
  \code{minnndist(X)} and \code{maxnndist(X)} are 
  equivalent to, but faster than, \code{min(nndist(X))}
  and \code{max(nndist(X))} respectively.

  The value is \code{NA} if \code{npoints(X) < 2}.
}
\value{
  A single numeric value (possibly \code{NA}).

  If \code{by} is given, the result is a numeric matrix
  giving the minimum or maximum nearest neighbour distance
  between each subset of \code{X}.
}
\seealso{
  \code{\link{nndist}}
}
\examples{
  min(nndist(swedishpines))
  minnndist(swedishpines)

  max(nndist(swedishpines))
  maxnndist(swedishpines)

  minnndist(lansing, positive=TRUE)

  if(interactive()) {
     X <- rpoispp(1e6)
     system.time(min(nndist(X)))
     system.time(minnndist(X))
  }

  minnndist(amacrine, by=marks(amacrine))
  maxnndist(amacrine, by=marks(amacrine))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
