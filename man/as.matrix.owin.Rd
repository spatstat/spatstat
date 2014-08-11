\name{as.matrix.owin}
\alias{as.matrix.owin}
\title{Convert Pixel Image to Matrix}
\description{
  Converts a pixel image to a matrix.
}
\usage{
  \method{as.matrix}{owin}(x, ...)
}
\arguments{
  \item{x}{A window (object of class \code{"owin"}).}
  \item{\dots}{Arguments passed to \code{\link{as.mask}} to control the
    pixel resolution.}
}
\details{
  The function \code{as.matrix.owin} converts a window to a
  logical matrux.

  It first converts the window \code{x} into a binary pixel mask
  using \code{\link{as.mask}}. It then extracts the pixel entries
  as a logical matrix.

  The resulting matrix has entries that are
  \code{TRUE} if the corresponding pixel is inside the window,
  and \code{FALSE} if it is outside.
  
  The function \code{as.matrix} is generic. The function
  \code{as.matrix.owin}
  is the method for windows (objects of class \code{"owin"}).

  Use \code{\link{as.im}} to convert a window to a pixel image.
}
\value{
  A logical matrix.
}
\examples{
  m <- as.matrix(letterR)
}
\seealso{
  \code{\link{as.matrix.im}},
  \code{\link{as.im}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
