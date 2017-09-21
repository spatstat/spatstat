\name{rotate.infline}
\alias{rotate.infline}
\alias{shift.infline}
\alias{reflect.infline}
\alias{flipxy.infline}
\title{
  Rotate or Shift Infinite Lines
}
\description{
  Given the coordinates of one or more infinite straight lines in the
  plane, apply a rotation or shift.
}
\usage{
\method{rotate}{infline}(X, angle = pi/2, \dots)

\method{shift}{infline}(X, vec = c(0,0), \dots)

\method{reflect}{infline}(X)

\method{flipxy}{infline}(X)
}
\arguments{
  \item{X}{
    Object of class \code{"infline"} representing one or more
    infinite straight lines in the plane.
}
  \item{angle}{
    Angle of rotation, in radians.
  }
  \item{vec}{
    Translation (shift) vector: a numeric vector of length 2,
    or a \code{list(x,y)}, or a point pattern containing one point.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  These functions are methods for the generic
  \code{\link{shift}}, \code{\link{rotate}},
  \code{\link{reflect}} and \code{\link{flipxy}}
  for the class \code{"infline"}.

  An object of class \code{"infline"} represents one or more
  infinite lines in the plane.
}
\value{
  Another \code{"infline"} object representing the
  result of the transformation.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{infline}}
}
\examples{
  L <- infline(v=0.5)

  plot(square(c(-1,1)), main="rotate lines", type="n")
  points(0, 0, pch=3)
  plot(L, col="green")
  plot(rotate(L, pi/12), col="red")
  plot(rotate(L, pi/6), col="red")
  plot(rotate(L, pi/4), col="red")

  L <- infline(p=c(0.4, 0.9), theta=pi* c(0.2, 0.6))

  plot(square(c(-1,1)), main="shift lines", type="n")
  L <- infline(p=c(0.7, 0.8), theta=pi* c(0.2, 0.6))
  plot(L, col="green")
  plot(shift(L, c(-0.5, -0.4)), col="red")

  plot(square(c(-1,1)), main="reflect lines", type="n")
  points(0, 0, pch=3)
  L <- infline(p=c(0.7, 0.8), theta=pi* c(0.2, 0.6))
  plot(L, col="green")
  plot(reflect(L), col="red")

}
\keyword{spatial}
\keyword{manip}
