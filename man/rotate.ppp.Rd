\name{rotate.ppp}
\alias{rotate.ppp}
\title{Rotate a Point Pattern}
\description{
  Rotates a point pattern
}
\usage{
 \method{rotate}{ppp}(X, angle=pi/2, \dots, centre=NULL)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
  \item{angle}{Angle of rotation.}
  \item{\dots}{
    Arguments passed to \code{\link{rotate.owin}} affecting
    the handling of the observation window, if it is a binary pixel mask.
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
  Another object of class \code{"ppp"} representing the
  rotated point pattern.
}
\details{
  The points of the pattern, and the window of observation, 
  are rotated about the origin by the angle specified.
  Angles are measured in
  radians, anticlockwise. The default is to rotate the pattern 90 degrees
  anticlockwise. If the points carry marks, these are preserved.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{rotate.owin}}
}
\examples{
  data(cells)
  X <- rotate(cells, pi/3)
  \dontrun{
  plot(X)
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
