\name{rotate.owin}
\alias{rotate.owin}
\title{Rotate a Window}
\description{
  Rotates a window
}
\usage{
 \method{rotate}{owin}(X, angle=pi/2, \dots, rescue=TRUE)
}
\arguments{
  \item{X}{A window (object of class \code{"owin"}).}
  \item{angle}{Angle of rotation.}
  \item{rescue}{
    Logical. If \code{TRUE}, the rotated window
    will be processed by \code{\link{rescue.rectangle}}.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}}
    controlling the resolution of the rotated window, if \code{X} is a
    binary pixel mask. Ignored if \code{X} is not a binary mask.
  }
}
\value{
  Another object of class \code{"owin"} representing the
  rotated window.
}
\details{
  Rotates the window by the specified angle. Angles are measured in
  radians, anticlockwise. The default is to rotate the window 90 degrees
  anticlockwise. The centre of rotation is the origin.
}
\seealso{
  \code{\link{owin.object}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  v <- rotate(w, pi/3)
  \dontrun{
  plot(v)
  }
  data(letterR)
  w <- as.mask(letterR)
  v <- rotate(w, pi/5)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
