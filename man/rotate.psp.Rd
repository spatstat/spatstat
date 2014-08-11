\name{rotate.psp}
\alias{rotate.psp}
\title{Rotate a Line Segment Pattern}
\description{
  Rotates a line segment pattern
}
\usage{
 \method{rotate}{psp}(X, angle=pi/2, \dots)
}
\arguments{
  \item{X}{A line segment pattern (object of class \code{"psp"}).}
  \item{angle}{Angle of rotation.}
  \item{\dots}{Arguments passed to \code{\link{rotate.owin}} affecting
    the handling of the observation window, if it is a binary pixel
    mask.
  }
}
\value{
  Another object of class \code{"psp"} representing the
  rotated line segment pattern.
}
\details{
  The line segments of the pattern, and the window of observation, 
  are rotated about the origin by the angle specified.
  Angles are measured in
  radians, anticlockwise. The default is to rotate the pattern 90 degrees
  anticlockwise. If the line segments carry marks, these are preserved.
}
\seealso{
  \code{\link{psp.object}},
  \code{\link{rotate.owin}},
  \code{\link{rotate.ppp}}
}
\examples{
  oldpar <- par(mfrow=c(2,1))
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  plot(X, main="original")
  Y <- rotate(X, pi/4)
  plot(Y, main="rotated")
  par(oldpar)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
