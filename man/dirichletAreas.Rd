\name{dirichletAreas}
\alias{dirichletAreas}
\title{
  Compute Areas of Tiles in Dirichlet Tessellation
}
\description{
  Calculates the area of each tile in the 
  Dirichlet-Voronoi tessellation of a point pattern.
}
\usage{
dirichletAreas(X)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
}
\details{
  This is an efficient algorithm to calculate the areas
  of the tiles in the Dirichlet-Voronoi tessellation.

  If the window of \code{X} is a binary pixel mask, the
  tile areas are computed by counting pixels. Otherwise the
  areas are computed exactly using analytic geometry.

  If any points of \code{X} are duplicated, the duplicates will
  have tile area zero.
}
\value{
  Numeric vector with one entry for each point of \code{X}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
\code{\link{dirichlet}},
\code{\link{dirichletVertices}}
}
\examples{
 aa <- dirichletAreas(cells)
}
\keyword{spatial}
\keyword{math}
\keyword{manip}
