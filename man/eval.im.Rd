\name{eval.im}
\alias{eval.im}
\title{Evaluate Expression Involving Pixel Images}
\description{
  Evaluates any expression involving one or more pixel images,
  and returns a pixel image.
}
\usage{
  eval.im(expr, envir, harmonize=TRUE, warn=TRUE)
}
\arguments{
  \item{expr}{An expression.}
  \item{envir}{Optional. The environment in which to evaluate the
    expression, or a named list containing pixel images to be used
    in the expression.}
  \item{harmonize}{
    Logical. Whether to resolve inconsistencies between
    the pixel grids.
  }
  \item{warn}{
    Logical. Whether to issue a warning if the pixel grids were
    inconsistent.
  }
}
\details{
  This function is a wrapper to make it easier to perform
  pixel-by-pixel calculations in an image. 

  Pixel images in \pkg{spatstat}
  are represented by objects of class \code{"im"}
  (see \code{\link{im.object}}). These are essentially matrices of
  pixel values, with extra attributes recording the pixel dimensions,
  etc.

  Suppose \code{X} is a pixel image. Then \code{eval.im(X+3)}
  will add 3 to the value of every pixel in \code{X}, and return
  the resulting pixel image.

  Suppose \code{X} and \code{Y} are two pixel images with compatible
  dimensions: they have the same number of pixels, the same physical
  size of pixels, and the same bounding box. Then
  \code{eval.im(X + Y)} will add the corresponding pixel values in
  \code{X} and \code{Y}, and return the resulting pixel image.

  In general, \code{expr} can be any expression in the R language involving
  (a) the \emph{names} of pixel images, (b) scalar
  constants, and (c) functions which are vectorised.
  See the Examples.

  First \code{eval.im} determines which of the \emph{variable names}
  in the expression \code{expr} refer to pixel images. Each such name
  is replaced by a matrix containing the pixel values. The expression is
  then evaluated. The result should be a matrix; it is taken as
  the matrix of pixel values.
  
  The expression \code{expr} must be vectorised.
  There must be at least one pixel image in the expression.

  All images must have compatible dimensions.
  If \code{harmonize=FALSE}, images that are incompatible will cause an error.
  If \code{harmonize=TRUE}, images that have incompatible dimensions
  will be resampled so that they are compatible; if \code{warn=TRUE},
  a warning will be issued.
}
\value{
  An image object of class \code{"im"}.
}
\seealso{
  \code{\link{im.apply}} for operations similar to \code{\link{apply}},
  such as taking the sum of a list of images.
  
  \code{\link{as.im}},
  \code{\link{compatible.im}},
  \code{\link{harmonise.im}},
  \code{\link{im.object}}
}
\examples{
  # test images
  X <- as.im(function(x,y) { x^2 - y^2 }, unit.square())
  Y <- as.im(function(x,y) { 3 * x + y }, unit.square())

  eval.im(X + 3)
  eval.im(X - Y)
  eval.im(abs(X - Y))
  Z <- eval.im(sin(X * pi) + Y)

  ## Use of 'envir'
  W <- eval.im(sin(U), list(U=density(cells)))
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
