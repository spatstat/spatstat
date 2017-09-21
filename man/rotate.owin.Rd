\name{rotate.owin}
\alias{rotate.owin}
\title{Rotate a Window}
\description{
  Rotates a window
}
\usage{
 \method{rotate}{owin}(X, angle=pi/2, \dots, rescue=TRUE, centre=NULL)
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
  \item{centre}{
    Centre of rotation.
    Either a vector of length 2, or a character string
    (partially matched to \code{"centroid"}, \code{"midpoint"}
    or \code{"bottomleft"}).
    The default is the coordinate origin \code{c(0,0)}.
  }
}
\value{
  Another object of class \code{"owin"} representing the
  rotated window.
}
\details{
  Rotates the window by the specified angle. Angles are measured in
  radians, anticlockwise. The default is to rotate the window 90 degrees
  anticlockwise. The centre of rotation is the origin, by default,
  unless \code{centre} is specified.
}
\seealso{
  \code{\link{owin.object}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  v <- rotate(w, pi/3)
  e <- rotate(w, pi/2, centre="midpoint")
  \dontrun{
  plot(v)
  }
  w <- as.mask(letterR)
  v <- rotate(w, pi/5)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
