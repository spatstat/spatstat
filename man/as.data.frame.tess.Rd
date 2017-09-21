\name{as.data.frame.tess}
\alias{as.data.frame.tess}
\title{Convert Tessellation to Data Frame}
\description{
  Converts a spatial tessellation object to a data frame.
}
\usage{
\method{as.data.frame}{tess}(x, \dots)
}
\arguments{
  \item{x}{
    Tessellation (object of class \code{"tess"}).
  }
  \item{\dots}{Further arguments passed to
    \code{\link{as.data.frame.owin}} or
    \code{\link{as.data.frame.im}} and ultimately to
    \code{\link[base:as.data.frame]{as.data.frame.default}}
    to determine the row names and other features.
  }
}
\details{
  This function converts the tessellation \code{x} to a data frame.

  If \code{x} is a pixel image tessellation (a pixel image with factor
  values specifying the tile membership of each pixel) then this
  pixel image is converted to a data frame by
  \code{\link{as.data.frame.im}}. The result is a data frame with
  columns \code{x} and \code{y} giving the pixel coordinates,
  and \code{Tile} identifying the tile containing the pixel.

  If \code{x} is a tessellation consisting of a rectangular grid of tiles
  or a list of polygonal tiles, then each tile is converted to a
  data frame by \code{\link{as.data.frame.owin}}, and these data frames
  are joined together, yielding a single large data frame containing
  columns \code{x}, \code{y} giving the coordinates of vertices of the
  polygons, and \code{Tile} identifying the tile.
}
\value{
  A data frame with columns named \code{x}, \code{y}, \code{Tile},
  and possibly other columns.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{as.data.frame.owin}},
  \code{\link{as.data.frame.im}}
}
\examples{
  Z <- as.data.frame(dirichlet(cells))
  head(Z, 10)
}
\keyword{spatial}
\keyword{methods}

