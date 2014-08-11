\name{gridweights}
\alias{gridweights}
\title{Compute Quadrature Weights Based on Grid Counts}
\description{
  Computes quadrature weights for a given set of points,
  using the ``counting weights'' for a grid of rectangular tiles.
}
\usage{
 gridweights(X, ntile, \dots, window=NULL, verbose=FALSE, npix=NULL, areas=NULL)
}
\arguments{
  \item{X}{Data defining a point pattern.}
  \item{ntile}{Number of tiles
    in each row and column of the rectangular grid.
    An integer vector of length 1 or 2.
  }
  \item{\dots}{Ignored.}
  \item{window}{Default window for the point pattern}
  \item{verbose}{Logical flag. If \code{TRUE}, information will be printed
    about the computation of the grid weights.
  }
  \item{npix}{Dimensions of pixel grid to use when
    computing a digital approximation to the tile areas.
  }
  \item{areas}{Vector of areas of the tiles, if they are already known.}
}
\value{
  Vector of nonnegative weights for each point in \code{X}.
}
\details{
  This function computes a set of quadrature weights
  for a given pattern of points
  (typically comprising both ``data'' and `dummy'' points).
  See \code{\link{quad.object}} for an explanation of quadrature
  weights and quadrature schemes.

  The weights are computed by the ``counting weights'' rule
  based on a regular grid of rectangular tiles.
  First \code{X} and (optionally) \code{window} are converted into a
  point pattern object. Then the bounding rectangle of the window of
  the point pattern is
  divided into a regular \code{ntile[1] * ntile[2]} grid of rectangular tiles.
  The weight attached to a point of \code{X} is the area of the tile
  in which it lies, divided by the number of points of \code{X} lying in
  that tile.

  For non-rectangular windows the tile areas are currently calculated
  by approximating the window as a binary mask. The accuracy of this
  approximation is controlled by \code{npix}, which becomes
  the argument \code{dimyx} of \code{\link{as.mask}}. 
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{dirichlet.weights}}
}
\examples{
  Q <- quadscheme(runifpoispp(10))
  X <- as.ppp(Q) # data and dummy points together
  w <- gridweights(X, 10)
  w <- gridweights(X, c(10, 10))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
