\name{tileindex}
\alias{tileindex}
\title{
  Determine Which Tile Contains Each Given Point
}
\description{
  Given a tessellation and a list of spatial points,
  determine which tile of the tessellation contains each of the
  given points.
}
\usage{
tileindex(x, y, Z)
}
\arguments{
  \item{x,y}{
    Spatial coordinates.
    Numeric vectors of equal length.
    (Alternatively \code{y} may be missing
    and \code{x} may be an object containing
    spatial coordinates).
  }
  \item{Z}{
    A tessellation (object of class \code{"tess"}).
  }
}
\details{
  This function determines which tile of the tessellation \code{Z}
  contains each of the spatial points
  with coordinates \code{(x[i],y[i])}.

  The result is a factor, of the same length as \code{x} and \code{y},
  indicating which tile contains each point. The levels of the factor
  are the names of the tiles of \code{Z}.
  Values are \code{NA} if the corresponding point lies outside the tessellation.
}
\value{
  A factor, of the same length as \code{x} and \code{y},
  whose levels are the names of the tiles of \code{Z}.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{cut.ppp}} and \code{\link{split.ppp}} to
  divide up the points of a point pattern according to
  a tessellation.

  \code{\link{as.function.tess}} to create a function whose
  value is the tile index.
}
\examples{
  X <- runifpoint(7)
  V <- dirichlet(X)
  tileindex(0.1, 0.4, V)
  tileindex(list(x=0.1, y=0.4), Z=V)
  tileindex(X, Z=V)
}
\keyword{spatial}
\keyword{manip}
