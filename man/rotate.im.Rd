\name{rotate.im}
\alias{rotate.im}
\title{Rotate a Pixel Image}
\description{
  Rotates a pixel image
}
\usage{
 \method{rotate}{im}(X, angle=pi/2, \dots)
}
\arguments{
  \item{X}{A pixel image (object of class \code{"im"}).}
  \item{angle}{Angle of rotation.}
  \item{\dots}{Ignored.}
}
\value{
  Another object of class \code{"im"} representing the
  rotated pixel image.
}
\details{
  The image is rotated about the origin by the angle specified.
  Angles are measured in radians, anticlockwise.
  The default is to rotate the image 90 degrees anticlockwise. 
}
\seealso{
  \code{\link{affine.im}},
  \code{\link{shift.im}},
  \code{\link{rotate}}
}
\examples{
  Z <- distmap(letterR)
  X <- rotate(Z)
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
