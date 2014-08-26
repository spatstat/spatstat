\name{raster.x}
\alias{raster.x}
\alias{raster.y}
\alias{raster.xy}
\title{Cartesian Coordinates for a Pixel Raster}
\description{
  Return the \eqn{x} and \eqn{y} coordinates
  of each pixel in a pixel image or binary mask.
}
\usage{
 raster.x(w, drop=FALSE)
 raster.y(w, drop=FALSE)
 raster.xy(w, drop=FALSE)
}
\arguments{
  \item{w}{
    A pixel image (object of class \code{"im"}) or a 
    mask window (object of class \code{"owin"}
    of type \code{"mask"}).
  }
  \item{drop}{
    Logical.
    If \code{TRUE}, then coordinates of pixels that lie outside the
    window are removed.
    If \code{FALSE} (the default) then the coordinates of every
    pixel in the containing rectangle are retained.
  }
}
\value{
  \code{raster.xy} returns a list with components \code{x} and \code{y}
  which are numeric vectors of equal length containing the pixel coordinates.

  If \code{drop=FALSE}, 
  \code{raster.x} and \code{raster.y} return 
  a matrix of the same dimensions as the pixel grid in \code{w},
  and giving the value of the \eqn{x} (or \eqn{y}) coordinate
  of each pixel in the raster.

  If \code{drop=TRUE}, 
  \code{raster.x} and \code{raster.y} return numeric vectors.
}
\details{
  The argument \code{w} should be either
  a pixel image (object of class \code{"im"})
  or a mask window (an object of class
  \code{"owin"} of type \code{"mask"}).

  If \code{drop=FALSE} (the default), the 
  functions \code{raster.x} and \code{raster.y} return
  a matrix of the same dimensions as the
  pixel image or mask itself, with entries giving the \eqn{x} coordinate
  (for \code{raster.x}) or \eqn{y} coordinate (for \code{raster.y})
  of each pixel in the pixel grid.

  If \code{drop=TRUE}, pixels that lie outside the
  window \code{w} (or outside the domain of the image \code{w})
  are removed, and \code{raster.x} and \code{raster.y} 
  return numeric vectors containing the coordinates of the
  pixels that are inside the window \code{w}.

  The function \code{raster.xy} returns a list
  with components \code{x} and \code{y}
  which are numeric vectors of equal length containing the pixel coordinates.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.mask}},
  \code{\link{pixelcentres}}
}
\examples{
  u <- owin(c(-1,1),c(-1,1)) # square of side 2
  w <- as.mask(u, eps=0.01) # 200 x 200 grid
  X <- raster.x(w)
  Y <- raster.y(w)
  disc <- owin(c(-1,1), c(-1,1), mask=(X^2 + Y^2 <= 1))
  \dontrun{plot(disc)}
  # approximation to the unit disc
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{manip}
