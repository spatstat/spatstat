\name{pixellate.psp}  %DontDeclareMethods
\alias{pixellate.psp}
\title{
  Convert Line Segment Pattern to Pixel Image
}
\description{
  Converts a line segment pattern to a pixel image
  by measuring the length of lines intersecting each pixel.
}
\usage{
\method{pixellate}{psp}(x, W=NULL, ..., weights = NULL)
}
\arguments{
  \item{x}{
    Line segment pattern (object of class \code{"psp"}).
  }
  \item{W}{
    Optional window (object of class \code{"owin"}) determining
    the pixel resolution.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution.
  }
  \item{weights}{
    Optional vector of weights associated with each line segment.
  }
}
\details{
  This function converts a line segment pattern to a pixel image by computing,
  for each pixel, the total length of
  intersection between the pixel and the line segments.
  
  This is a method for the generic function
  \code{\link{pixellate}} for the class of line segment patterns.

  The pixel raster is determined by \code{W}
  and the optional arguments \code{\dots}.
  If \code{W} is missing or \code{NULL}, it defaults to the window
  containing \code{x}.
  Then \code{W} is converted to a
  binary pixel mask using \code{\link{as.mask}}. The arguments
  \code{\dots} are passed to \code{\link{as.mask}} to
  control the pixel resolution.

  If \code{weights} are given, then the length of the intersection
  between line segment \code{i} and pixel \code{j} is multiplied by
  \code{weights[i]} before the lengths are summed for each pixel.
}
\value{
  A pixel image (object of class \code{"im"}) with numeric values.
}
\seealso{
  \code{\link{pixellate}},
  \code{\link{as.mask}},
  \code{\link{as.mask.psp}}.
  
  Use \code{\link{as.mask.psp}} if you only want to know
  which pixels are intersected by lines.
}
\examples{
  X <- psp(runif(10),runif(10), runif(10), runif(10), window=owin())
  plot(pixellate(X))
  plot(X, add=TRUE)
  sum(lengths.psp(X))
  sum(pixellate(X))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}

