\name{which.max.im}
\alias{which.max.im}
\title{Identify Pixelwise Maximum of Several Pixel Images}
\description{
  Given a list of pixel images, this function 
  identifies the image that has the largest value
  at each pixel, and returns an image.
}
\usage{
  which.max.im(x)
}
\arguments{
  \item{x}{A list of images (objects of class \code{"im"}).}
}
\details{
  \code{x} should be a list of pixel images.
  All images must have compatible dimensions.

  For each pixel, the algorithm identifies which of the
  images in the list \code{x} has the largest value at that pixel.
  The index of this image becomes the pixel value in the output image.
  If \code{names(x)} is not null, then the indices are replaced by
  these names.
}
\value{
  An image (object of class \code{"im"}) with factor values.
}
\seealso{
  \code{\link{eval.im}},
  \code{\link{im.object}}
}
\examples{
  # test images
  X <- as.im(function(x,y) { x^2 - y^2 }, unit.square())
  Y <- as.im(function(x,y) { x - y }, unit.square())
  which.max.im(list(X=X,Y=Y))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
