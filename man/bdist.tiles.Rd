\name{bdist.tiles}
\alias{bdist.tiles}
\title{Distance to Boundary of Window}
\description{
  Computes the shortest distances 
  from each tile in a tessellation
  to the boundary of the window.
}
\usage{
 bdist.tiles(X)
}
\arguments{
  \item{X}{A tessellation (object of class \code{"tess"}).}
}
\value{
  A numeric vector, 
  giving the shortest distance from each tile in the tessellation
  to the boundary of the window.
  Entries of the vector correspond to the entries of \code{tiles(X)}.
}
\details{
  This function computes, for each tile \eqn{s_i}{s[[i]]}
  in the tessellation \code{X}, the shortest distance
  from \eqn{s_i}{s[[i]]}
  to the boundary of the window \eqn{W} containing the tessellation.
}
\seealso{
  \code{\link{tess}},
  \code{\link{bdist.points}},
  \code{\link{bdist.pixels}}
}
\examples{
  P <- runifpoint(15)
  X <- dirichlet(P)
  plot(X, col="red")
  B <- bdist.tiles(X)
  # identify tiles that do not touch the boundary
  plot(X[B > 0], add=TRUE, col="green", lwd=3)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
