\name{dirichlet}
\alias{dirichlet}
\title{Dirichlet Tessellation of Point Pattern}
\description{
  Computes the Dirichlet tessellation of a spatial point pattern.
  Also known as the Voronoi or Thiessen tessellation.
}
\usage{
dirichlet(X)
}
\arguments{
  \item{X}{Spatial point pattern (object of class \code{"ppp"}).}
}
\details{
  In a spatial point pattern \code{X}, the Dirichlet tile associated
  with a particular point \code{X[i]} is the region of space that is
  closer to \code{X[i]} than to any other point in \code{X}. The
  Dirichlet tiles divide the two-dimensional plane into disjoint
  regions, forming a tessellation.

  The Dirichlet tessellation is also known as the Voronoi or
  Thiessen tessellation.

  This function computes the Dirichlet tessellation (within the original
  window of \code{X}) using the function \code{\link[deldir]{deldir}}
  in the package \pkg{deldir}.

  To ensure that there is a one-to-one correspondence between the
  points of \code{X} and the tiles of \code{dirichlet(X)},
  duplicated points in \code{X} should first be removed by
  \code{X <- unique(X, rule="deldir")}. 
  
  The tiles of the tessellation will be computed as polygons
  if the original window is a rectangle or a polygon.
  Otherwise the tiles will be computed as binary masks. 
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\seealso{
  \code{\link{tess}},
  \code{\link{delaunay}},
  \code{\link{ppp}}
}
\examples{
  X <- runifpoint(42)
  plot(dirichlet(X))
  plot(X, add=TRUE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
