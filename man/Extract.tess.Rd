\name{Extract.tess}
\alias{[.tess}
\alias{[<-.tess}
\title{Extract or Replace Subset of Tessellation}
\description{
  Extract, change or delete a subset of the tiles of a tessellation,
  to make a new tessellation.
}
\usage{
  \method{[}{tess}(x, i, \dots)
  \method{[}{tess}(x, i, \dots) <- value
}
\arguments{
  \item{x}{A tessellation (object of class \code{"tess"}).}
  \item{i}{
     Subset index for the tiles of the tessellation.
     Alternatively a window (object of class \code{"owin"}).
  }
  \item{\dots}{
     One argument that specifies the subset to be extracted or changed.
     Any valid format for the subset index in a list.
  }
  \item{value}{
     Replacement value for the selected tiles of the tessellation.
     A list of windows (objects of class \code{"owin"}) or \code{NULL}. 
  }
}
\details{
  A tessellation (object of class \code{"tess"}, see \code{\link{tess}})
  is effectively a list of tiles (spatial regions) that cover a spatial region.
  The subset operator \code{[.tess} extracts some
  of these tiles and forms a new tessellation, which of course covers a
  smaller region than the original. 

  For \code{[.tess} only, the subset index can also be a window
  (object of class \code{"owin"}). The tessellation \code{x}
  is then intersected with the window.

  The replacement operator changes the selected tiles. The replacement
  \code{value} may be either \code{NULL} (which causes the selected tiles
  to be removed from \code{x}) or a list of the same length as
  the selected subset. The entries of \code{value} may be windows
  (objects of class \code{"owin"}) or \code{NULL} to indicate that the
  corresponding tile should be deleted.

  Generally it does not make sense to replace a tile in a tessellation
  with a completely different tile, because the tiles are expected to
  fit together. However this facility is sometimes useful for making
  small adjustments to polygonal tiles.
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\seealso{
  \code{\link{tess}}, \code{\link{tiles}}, \code{\link{intersect.tess}}.
}
\examples{
   \testonly{op <- spatstat.options(npixel=10)}
   A <- tess(xgrid=0:4, ygrid=0:3)
   B <- A[c(1, 3, 7)]
   E <- A[-1]
   A[c(2, 5, 11)] <- NULL
   \testonly{spatstat.options(op)}
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{manip}
