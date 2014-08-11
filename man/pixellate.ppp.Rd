\name{pixellate.ppp}
\alias{pixellate.ppp}
\alias{as.im.ppp}
\title{Convert Point Pattern to Pixel Image}
\description{
  Converts a point pattern to a pixel image. The value in each pixel
  is the number of points falling in that pixel, and is typically
  either 0 or 1.
}
\usage{
\method{pixellate}{ppp}(x, W=NULL, \dots, weights = NULL, padzero=FALSE)

\method{as.im}{ppp}(X, \dots)
}
\arguments{
  \item{x,X}{Point pattern (object of class \code{"ppp"}).}
  \item{\dots}{Arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution}
  \item{W}{Optional window mask (object of class \code{"owin"}) determining
    the pixel raster.
  }
  \item{weights}{Optional vector of weights associated with the points.}
  \item{padzero}{Logical flag indicating whether to set pixel values
    to zero outside the window. 
  }
}
\details{
  The functions \code{pixellate.ppp} and \code{as.im.ppp}
  convert a spatial point pattern \code{x} into a pixel
  image, by counting the number of points (or the total weight of
  points) falling in each pixel. 

  Calling \code{as.im.ppp} is equivalent to
  calling \code{pixellate.ppp} with its default arguments.
  Note that \code{pixellate.ppp} is more general than \code{as.im.ppp}
  (it has additional arguments for greater flexibility).

  The functions \code{as.im.ppp} and \code{pixellate.ppp}
  are methods for the generic functions \code{\link{as.im}}
  and \code{\link{pixellate}} respectively,
  for the class of point patterns.

  The pixel raster (in which points are counted) is determined
  by the argument \code{W} if it is present (for \code{pixellate.ppp} only).
  In this case \code{W} should be a binary mask (a window object of
  class \code{"owin"} with type \code{"mask"}).
  Otherwise the pixel raster is determined by
  extracting the window containing \code{x} and converting it to a
  binary pixel mask using \code{\link{as.mask}}. The arguments
  \code{\dots} are passed to \code{\link{as.mask}} to
  control the pixel resolution.

  If \code{weights} is \code{NULL}, then for each pixel
  in the mask, the algorithm counts how many points in \code{x} fall
  in the pixel. This count is usually either 0 (for a pixel with no data
  points in it) or 1 (for a pixel containing one data point) but may be
  greater than 1. The result is an image with these counts as its pixel values.

  If \code{weights} is given, it should be a numeric vector of the same
  length as the number of points in \code{x}. For each pixel, the
  algorithm finds the total weight associated with points in \code{x} that fall
  in the given pixel. The result is an image with these total weights
  as its pixel values.

  By default (if \code{zeropad=FALSE}) the resulting pixel image has the same
  spatial domain as the window of the point pattern \code{x}. If
  \code{zeropad=TRUE} then the resulting pixel image has a rectangular
  domain; pixels outside the original window are assigned the value zero.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\seealso{
  \code{\link{pixellate}},
  \code{\link{im}}, \code{\link{as.im}},
  \code{\link{density.ppp}}, \code{\link{smooth.ppp}}.
}
\examples{
  data(humberside)
  plot(pixellate(humberside))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
