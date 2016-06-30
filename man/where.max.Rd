\name{where.max}
\alias{where.max}
\alias{where.min}
\title{
  Find Location of Maximum in a Pixel Image
}
\description{
  Finds the spatial location(s) where a given pixel image
  attains its maximum or minimum value.
}
\usage{
  where.max(x, first = TRUE)
  where.min(x, first = TRUE)
}
\arguments{
  \item{x}{
    A pixel image (object of class \code{"im"}).
  }
  \item{first}{
    Logical value. If \code{TRUE} (the default), then only one location
    will be returned. If \code{FALSE}, then all locations where the
    maximum is achieved will be returned.
  }
}
\details{
  This function finds the spatial location or locations where the
  pixel image \code{x} attains its maximum or minimum value.
  The result is a point pattern giving the locations.

  If \code{first=TRUE} (the default), then only one location will
  be returned, namely the location with the smallest \eqn{y} coordinate
  value which attains the maximum or minimum.
  This behaviour is analogous to the functions
  \code{\link[base]{which.min}} and
  \code{\link[base:which.min]{which.max}}.

  If \code{first=FALSE}, then the function returns
  the locations of all pixels where the
  maximum (or minimum) value is attained. This could be a large
  number of points.
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{max.im}},
  \code{\link{eval.im}},
  \code{\link{solutionset}}
}
\examples{
   D <- distmap(letterR, invert=TRUE)
   plot(D)
   plot(where.max(D), add=TRUE, pch=16, cols="green")
}
\keyword{spatial}
\keyword{math}
