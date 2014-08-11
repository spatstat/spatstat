\name{pixelquad}
\alias{pixelquad}
\title{Quadrature Scheme Based on Pixel Grid}
\description{
  Makes a quadrature scheme with a dummy point at every pixel
  of a pixel image.
}
\usage{
pixelquad(X, W = as.owin(X))
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"}) containing the
    data points for the quadrature scheme.
  }
  \item{W}{
    Specifies the pixel grid.
    A pixel image (object of class \code{"im"}),
    a window (object of class \code{"owin"}), or anything that can
    be converted to a window by \code{\link{as.owin}}.
  }
}
\value{
  An object of class \code{"quad"} describing the quadrature scheme
  (data points, dummy points, and quadrature weights)
  suitable as the argument \code{Q} of the function \code{\link{ppm}()} for
  fitting a point process model.

  The quadrature scheme can be inspected using the
  \code{print} and \code{plot} methods for objects
  of class \code{"quad"}. 
}
\details{
  This is a method for producing a quadrature scheme
  for use by \code{\link{ppm}}. It is an alternative to
  \code{\link{quadscheme}}.
  
  The function \code{\link{ppm}} fits a point process model to an
  observed point pattern using 
  the Berman-Turner quadrature approximation (Berman and Turner, 1992;
  Baddeley and Turner, 2000) to the pseudolikelihood of the model. 
  It requires a quadrature scheme consisting of 
  the original data point pattern, an additional pattern of dummy points,
  and a vector of quadrature weights for all these points.
  Such quadrature schemes are represented by objects of class
  \code{"quad"}. See \code{\link{quad.object}} for a description of this class.

  Given a grid of pixels, this function creates a quadrature scheme
  in which there is one dummy point at the centre of each pixel. The
  counting weights are used (the weight attached to each quadrature
  point is 1 divided by the number of quadrature points falling
  in the same pixel).

  The argument \code{X} specifies the locations of the data points
  for the quadrature scheme. Typically this would be a point pattern
  dataset.

  The argument \code{W} specifies the grid of pixels for the dummy
  points of the quadrature scheme. It should be a pixel image
  (object of class \code{"im"}), a window (object of class
  \code{"owin"}), or anything that can
  be converted to a window by \code{\link{as.owin}}. If \code{W} is a
  pixel image or a binary mask (a window of type \code{"mask"})
  then the pixel grid of \code{W} will be used. If \code{W} is a
  rectangular or polygonal window, then it will first be converted to a
  binary mask using \code{\link{as.mask}} at the default pixel
  resolution.
}
\examples{
  W <- owin(c(0,1),c(0,1))
  X <- runifpoint(42, W)
  W <- as.mask(W,dimyx=128)
  pixelquad(X,W)
}
\seealso{
  \code{\link{quadscheme}},
  \code{\link{quad.object}},
  \code{\link{ppm}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
