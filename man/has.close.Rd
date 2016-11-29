\name{has.close}
\alias{has.close}
\alias{has.close.default}
\alias{has.close.ppp}
\alias{has.close.pp3}
\title{
  Check Whether Points Have Close Neighbours
}
\description{
  For each point in a point pattern, determine whether the
  point has a close neighbour in the same pattern.
}
\usage{
  has.close(X, r, Y=NULL, \dots)

  \method{has.close}{default}(X,r, Y=NULL, \dots, periodic=FALSE)

  \method{has.close}{ppp}(X,r, Y=NULL, \dots, periodic=FALSE, sorted=FALSE)

  \method{has.close}{pp3}(X,r, Y=NULL, \dots, periodic=FALSE, sorted=FALSE)
}
\arguments{
  \item{X,Y}{
    Point patterns of class \code{"ppp"} or \code{"pp3"} or \code{"lpp"}.
  }
  \item{r}{
    Threshold distance: a number greater than zero.
  }
  \item{periodic}{
    Logical value indicating whether to measure distances in the
    periodic sense, so that opposite sides of the (rectangular) window
    are treated as identical.
  }
  \item{sorted}{
    Logical value, indicating whether the points of \code{X}
    (and \code{Y}, if given) are already sorted into increasing order of the
    \eqn{x} coordinates.
  }
  \item{\dots}{Other arguments are ignored.}
}
\details{
  This is simply a faster version of \code{(nndist(X) <= r)}
  or \code{(nncross(X,Y,what="dist") <= r)}.

  \code{has.close(X,r)} determines, for each point in the pattern \code{X},
  whether or not this point has a neighbour in the same pattern \code{X}
  which lies at a distance less than or equal to \code{r}.

  \code{has.close(X,r,Y)} determines, for each point in the pattern \code{X},
  whether or not this point has a neighbour in the \emph{other} pattern
  \code{Y} which lies at a distance less than or equal to \code{r}.

  The function \code{has.close} is generic, with methods for
  \code{"ppp"} and \code{"pp3"} and a default method.
}
\value{
  A logical vector, with one entry for each point of \code{X}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{nndist}}
}
\examples{
  has.close(redwood, 0.05)
  with(split(amacrine), has.close(on, 0.05, off))
}
\keyword{spatial}
\keyword{math}

