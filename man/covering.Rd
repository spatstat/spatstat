\name{covering}
\alias{covering}
\title{Cover Region with Discs}
\description{
  Given a spatial region, this function finds an efficient
  covering of the region using discs of a chosen radius.
}
\usage{
 covering(W, r, \dots, giveup=1000)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}).
  }
  \item{r}{positive number: the radius of the covering discs.}
  \item{\dots}{
    extra arguments passed to \code{\link{as.mask}}
    controlling the pixel resolution for the calculations.
  }
  \item{giveup}{
    Maximum number of attempts to place additional discs.
  }
}
\value{
  A point pattern (object of class \code{"ppp"})
  giving the centres of the discs.
}
\details{
  This function finds an efficient covering of the
  window \code{W} using discs of the given radius \code{r}.
  The result is a point pattern giving the centres of the discs.

  The algorithm tries to use as few discs as possible,
  but is not guaranteed to find the minimal number of discs.
  It begins by placing a hexagonal grid of points inside \code{W},
  then adds further points until every location inside \code{W} lies
  no more than \code{r} units away from one of the points.
}
\examples{
  rr <- 0.5
  X <- covering(letterR, rr)
  plot(grow.rectangle(Frame(X), rr), type="n", main="")
  plot(X, pch=16, add=TRUE, col="red")
  plot(letterR, add=TRUE, lwd=3)
  plot(X \%mark\% (2*rr), add=TRUE, markscale=1)
}
\author{
  \adrian
}
\keyword{spatial}
\keyword{math}
