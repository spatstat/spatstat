\name{harmonise.im}
\alias{harmonise.im}
\alias{harmonize.im}
\title{Make Pixel Images Compatible}
\description{
  Convert several pixel images to a common pixel raster.
}
\usage{
\method{harmonise}{im}(\dots)

\method{harmonize}{im}(\dots)
}
\arguments{
  \item{\dots}{
    Any number of pixel images (objects of class \code{"im"})
    or data which can be converted to pixel images by \code{\link{as.im}}.
  } 
}
\details{
  This function makes any number of pixel images compatible,
  by converting them all to a common pixel grid.

  The command \code{\link{harmonise}} is generic. This is the
  method for objects of class \code{"im"}.
  
  At least one of the arguments \code{\dots} must be a pixel image.
  Some arguments may be windows (objects of class \code{"owin"}),
  functions (\code{function(x,y)}) or numerical constants. These will be
  converted to images using \code{\link{as.im}}.
  
  The common pixel grid is determined by inspecting all the pixel
  images in the argument list, computing the bounding box of all the
  images, then finding the image with the highest spatial resolution, 
  and extending its pixel grid to cover the bounding box. 

  The return value is a list with entries corresponding to the input
  arguments.
  If the arguments were named (\code{name=value}) then the return value
  also carries these names.

  If you just want to determine the appropriate pixel resolution,
  without converting the images, use \code{\link{commonGrid}}.
}
\value{
  A list,
  of length equal to the number of arguments \code{\dots},
  whose entries are pixel images.
}
\author{\adrian
  
  
  and \rolf
  
}
\examples{
   A <- setcov(square(1))
   B <- function(x,y) { x }
   G <- density(runifpoint(42))
   harmonise(X=A, Y=B, Z=G)
}
\seealso{
  \code{\link{commonGrid}},
  \code{\link{compatible.im}},
  \code{\link{as.im}}
}
\keyword{spatial}
\keyword{manip}
