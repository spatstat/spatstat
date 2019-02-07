\name{lineartileindex}
\alias{lineartileindex}
\title{
  Determine Which Tile Contains Each Given Point on a Linear Network
}
\description{
  Given a tessellation on a linear network,
  and a list of points on the network,
  determine which tile of the tessellation contains each of the given points.
}
\usage{
lineartileindex(seg, tp, Z, method = c("encode", "C", "interpreted"))
}
\arguments{
  \item{seg,tp}{
    Vectors of local coordinates of the query points.
    See Details.
  }
  \item{Z}{
    A tessellation on a linear network (object of class
    \code{"lintess"}).
  }
  \item{method}{
    Internal use only.
  }
}
\details{
  This low-level function is 
  the analogue of \code{\link{tileindex}} for linear networks.
  For a tessellation \code{Z} on a linear network,
  and a list of query points on the same network,
  the function determines which tile of the tessellation 
  contains each query point.

  Argument \code{Z} should be a tessellation on a linear network
  (object of class \code{"lintess"}).

  The vectors \code{seg} and \code{tp} specify the locations
  of the query points, on the same network, using local coordinates:
  \code{seg} contains integer values specifying which
  segment of the network contains each query point;
  \code{tp} contains numeric values between 0 and 1 specifying the
  fractional position along that segment.

  The result is a factor, of the same length as \code{seg} and \code{tp},
  indicating which tile contains each point. The levels of the factor
  are the names of the tiles of \code{Z}.
}
\value{
  A factor, of the same length as \code{seg} and \code{tp},
  whose levels are the names of the tiles of \code{Z}.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{lintess}}.
  
  \code{\link{as.linfun.lintess}} to create a function whose
  value is the tile index.

  \code{\link{cut.lpp}} for a neater way to classify
  the points of a point pattern on a linear network
  according to a tessellation on the network.
}
\examples{
   Z <- lineardirichlet(runiflpp(15, simplenet))
   X <- runiflpp(10, simplenet)
   coX <- coords(X)
   ii <- lineartileindex(coX$seg, coX$tp, Z)
}
\keyword{spatial}
\keyword{manip}
