\name{spokes}
\alias{spokes}
\title{Spokes pattern of dummy points}
\description{
  Generates a pattern of dummy points in a window,
  given a data point pattern.
  The dummy points lie on the radii of circles emanating from each
  data point.
}
\usage{
 spokes(x, y, nrad = 3, nper = 3, fctr = 1.5, Mdefault = 1)
}
\arguments{
  \item{x}{
    Vector of \eqn{x} coordinates of data points, or a list
    with components \code{x} and \code{y}, or a point pattern
    (an object of class \code{ppp}).
  }
  \item{y}{
    Vector of \eqn{y} coordinates of data points.  Ignored
    unless \code{x} is a vector.
  }
  \item{nrad}{
    Number of radii emanating from each data point.
  }
  \item{nper}{
    Number of dummy points per radius.
  }
  \item{fctr}{
    Scale factor.
    Length of largest spoke radius is \code{fctr * M}
    where \code{M} is the mean nearest neighbour distance
    for the data points.
  }
  \item{Mdefault}{
    Value of \code{M} to be used if \code{x} has length 1.
  }
}
\value{
  If argument \code{x} is a point pattern, a point pattern with
  window equal to that of \code{x}.  Otherwise a list with two
  components \code{x} and \code{y}.  In either case the components
  \code{x} and \code{y} of the value are numeric vectors giving
  the coordinates of the dummy points.
}
\details{
  This function is useful in creating dummy points for quadrature
  schemes (see \code{\link{quadscheme}}).

  Given the data points, the function creates a collection of
  \code{nrad * nper * length(x)} dummy points.

  Around each data point \code{(x[i],y[i])} there are
  \code{nrad * nper} dummy points, lying on \code{nrad} radii
  emanating from \code{(x[i],y[i])}, with \code{nper} dummy points
  equally spaced along each radius.

  The (equal) spacing of dummy points along each radius is
  controlled by the factor \code{fctr}. 
  The distance from a data point to the furthest of its associated
  dummy points is \code{fctr * M}
  where \code{M} is the mean nearest neighbour distance
  for the data points.

  If there is only one data point the nearest neighbour distance
  is infinite, so the value \code{Mdefault} will be used in place
  of \code{M}.

  If \code{x} is a point pattern, then the value returned is
  also a point pattern, which is clipped to the window
  of \code{x}.  Hence there may be fewer than
  \code{nrad * nper * length(x)} dummy points in the pattern
  returned.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{quadscheme}},
  \code{\link{inside.owin}},
  \code{\link{gridcentres}},
  \code{\link{stratrand}}
}
\examples{
  dat <- runifrect(10)
  dum <- spokes(dat$x, dat$y, 5, 3, 0.7)
  plot(dum)
  Q <- quadscheme(dat, dum, method="dirichlet")
  plot(Q, tiles=TRUE)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{datagen}
