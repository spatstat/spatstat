\name{setcov}
\alias{setcov}
\title{Set Covariance of a Window}
\description{
  Computes the set covariance function of a window.
}
\usage{
 setcov(W, V=W, \dots)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}.
  }
  \item{V}{
    Optional. Another window.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}}
    to control the pixel resolution.
  }
}
\value{
  A pixel image (an object of class \code{"im"}) representing the
  set covariance function of \code{W},
  or the cross-covariance of \code{W} and \code{V}.
}
\details{
  The set covariance function of a region \eqn{W} in the plane
  is the function \eqn{C(v)} defined for each vector \eqn{v}
  as the area of the intersection between \eqn{W} and \eqn{W+v},
  where \eqn{W+v} is the set obtained by shifting (translating)
  \eqn{W} by \eqn{v}.

  We may interpret \eqn{C(v)} as the area of the set of
  all points \eqn{x} in \eqn{W} such that \eqn{x+v} also lies in
  \eqn{W}.
  
  This command computes a discretised approximation to
  the set covariance function of any
  plane region \eqn{W} represented as a window object (of class
  \code{"owin"}, see \code{\link{owin.object}}). The return value is
  a pixel image (object of class \code{"im"}) whose greyscale values
  are values of the set covariance function.

  The set covariance is computed using the Fast Fourier Transform,
  unless \code{W} is a rectangle, when an exact formula is used.

  If the argument \code{V} is present, then \code{setcov(W,V)}
  computes the set \emph{cross-covariance} function \eqn{C(x)}
  defined for each vector \eqn{x}
  as the area of the intersection between \eqn{W} and \eqn{V+x}.
}
\seealso{
  \code{\link{imcov}},
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{erosion}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  v <- setcov(w)
  plot(v)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
