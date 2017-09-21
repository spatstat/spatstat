\name{bdist.pixels}
\alias{bdist.pixels}
\title{Distance to Boundary of Window}
\description{
  Computes the distances 
  from each pixel in a window to the boundary of the window.
}
\usage{
 bdist.pixels(w, \dots, style="image", method=c("C", "interpreted"))
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}).}
  \item{\dots}{Arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution.}
  \item{style}{Character string determining the format of
    the output: either \code{"matrix"}, \code{"coords"} or
    \code{"image"}.
  }
  \item{method}{Choice of algorithm to use when \code{w} is polygonal.}
}
\value{
  If \code{style="image"}, a pixel image (object of class \code{"im"})
  containing the distances from each pixel in the image raster
  to the boundary of the window.

  If \code{style="matrix"},
  a matrix giving the distances from each pixel in the image raster
  to the boundary of the window. Rows of this matrix correspond to
  the \eqn{y} coordinate and columns to the \eqn{x} coordinate.
  
  If \code{style="coords"}, a list with three components
  \code{x,y,z}, where \code{x,y} are vectors of length \eqn{m,n}
  giving the \eqn{x} and \eqn{y} coordinates respectively,
  and \code{z} is an \eqn{m \times n}{m x n} matrix such that
  \code{z[i,j]} is the distance from \code{(x[i],y[j])} to the
  boundary of the window. Rows of this matrix correspond to the
  \eqn{x} coordinate and columns to the \eqn{y} coordinate.
  This result can be plotted with \code{persp}, \code{image}
  or \code{contour}.
}
\details{
  This function computes, for each pixel \eqn{u}
  in the window \code{w}, the shortest distance
  \eqn{d(u, W^c)}{dist(u, W')} from \eqn{u}
  to the boundary of \eqn{W}.

  If the window is a binary mask then the distance from each pixel
  to the boundary is computed using the distance transform algorithm
  \code{\link{distmap.owin}}. The result is equivalent to
  \code{distmap(W, invert=TRUE)}.

  If the window is a rectangle or a polygonal region,
  the grid of pixels is determined by the arguments \code{"\dots"} 
  passed to \code{\link{as.mask}}. The distance from each pixel to the
  boundary is calculated exactly, using analytic geometry.
  This is slower but more accurate than in the case of a binary mask.

  For software testing purposes, there are two implementations
  available when \code{w} is a polygon: the default is \code{method="C"}
  which is much faster than \code{method="interpreted"}.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{erosion}},
  \code{\link{bdist.points}},
  \code{\link{bdist.tiles}},
  \code{\link{distmap.owin}}.
}
\examples{
  u <- owin(c(0,1),c(0,1))
  d <- bdist.pixels(u, eps=0.01)
  image(d)
  d <- bdist.pixels(u, eps=0.01, style="matrix")
  mean(d >= 0.1)
  # value is approx (1 - 2 * 0.1)^2 = 0.64
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
