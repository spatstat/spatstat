\name{bdist.points}
\alias{bdist.points}
\title{Distance to Boundary of Window}
\description{
  Computes the distances 
  from each point of a point pattern
  to the boundary of the window.
}
\usage{
 bdist.points(X)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
}
\value{
  A numeric vector, giving the distances from each point of the pattern
  to the boundary of the window.
}
\details{
  This function computes, for each point \eqn{x_i}{x[i]}
  in the point pattern \code{X}, the shortest distance
  \eqn{d(x_i, W^c)}{dist(x[i], W')} from \eqn{x_i}{x[i]}
  to the boundary of the window \eqn{W} of observation.

  If the window \code{Window(X)} is of type \code{"rectangle"}
  or \code{"polygonal"}, then these distances are computed by
  analytic geometry and are exact, up to rounding errors.
  If the window is of type \code{"mask"} then the distances
  are computed using the real-valued distance transform,
  which is an approximation with maximum error equal to the width
  of one pixel in the mask.
}
\seealso{
  \code{\link{bdist.pixels}},
  \code{\link{bdist.tiles}},
  \code{\link{ppp.object}},
  \code{\link{erosion}}
}
\examples{
  data(cells)
  d <- bdist.points(cells)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
