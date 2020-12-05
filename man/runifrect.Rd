\name{runifrect}
\alias{runifrect}
\title{Generate N Uniform Random Points in a Rectangle}
\description{
  Generate a random point pattern,
  containing \eqn{n} independent uniform random points,
  inside a specified rectangle.
}
\usage{
runifrect(n, win = owin(c(0, 1), c(0, 1)), nsim = 1, drop = TRUE)
}
\arguments{
  \item{n}{
    Number of points.
  }
  \item{win}{
    Rectangular window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}},
    which must specify a rectangle.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  This function is a slightly faster version of
  \code{\link[spatstat]{runifpoint}}
  for the special case where the window is a rectangle.

  The function generates \code{n} independent random points,
  uniformly distributed in the window \code{win},
  by assigning uniform random values to the cartesian coordinates.

  For normal usage we recommend \code{\link[spatstat]{runifpoint}}
  because it is more flexible.
  However, \code{runifrect} is slightly faster (when the window is a rectangle),
  and may be preferable in very computationally-demanding tasks.
}
\value{
  A point pattern (an object of class \code{"ppp"})
  if \code{nsim=1} and \code{drop=TRUE},
  otherwise a list of point patterns.
}
\seealso{
\code{\link{ppp.object}},
\code{\link{owin.object}},
\code{\link[spatstat]{runifpoint}},
\code{\link[spatstat]{rpoispp}},
\code{\link[spatstat]{rpoint}}
}
\examples{
 # 42 random points in the unit square
 pp <- runifrect(42)
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{datagen}

