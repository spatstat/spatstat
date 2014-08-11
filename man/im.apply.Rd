\name{im.apply}
\alias{im.apply}
\title{
  Apply Function Pixelwise to List of Images 
}
\description{
  Returns a pixel image obtained by applying a function
  to the values of corresponding pixels in several pixel images.
}
\usage{
im.apply(X, FUN, ...)
}
\arguments{
  \item{X}{
    A list of pixel images (objects of class \code{"im"}).
  }
  \item{FUN}{
    A function that can be applied to vectors,
    or a character string giving the name of such a function.
  }
  \item{\dots}{
    Additional arguments to \code{FUN}.
  }
}
\details{
  The argument \code{X} should be a list of pixel images
  (objects of class \code{"im"}). If the images do not have
  identical pixel grids, they will be converted to a common
  grid using \code{\link{harmonise.im}}.

  At each pixel location, the values of the images in \code{X}
  at that pixel will be extracted as a vector. The function
  \code{FUN} will be applied to this vector. The result (which should be
  a single value) becomes the pixel value of the resulting image.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\seealso{
  \code{\link{eval.im}} for algebraic operations with images.
}
\examples{
  DA <- density(split(amacrine))
  DA
  im.apply(DA, max)
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
