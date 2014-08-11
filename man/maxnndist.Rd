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
minnndist(X)
maxnndist(X)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
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
}
\seealso{
  \code{\link{nndist}}
}
\examples{
  min(nndist(swedishpines))
  minnndist(swedishpines)

  max(nndist(swedishpines))
  maxnndist(swedishpines)

  if(interactive()) {
     X <- rpoispp(1e6)
     system.time(min(nndist(X)))
     system.time(minnndist(X))
  }
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner and Ege Rubak.
}
\keyword{spatial}
\keyword{math}
