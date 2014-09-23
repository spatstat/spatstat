\name{dirichlet.vertices}
\alias{dirichlet.vertices}
\title{
   Vertices of Dirichlet Tessellation
}
\description{
  Computes the Dirichlet-Voronoi tessellation of a point pattern
  and extracts the vertices of the tiles.
}
\usage{
dirichlet.vertices(X)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
}
\details{
  The Dirichlet tessellation is described in the help for

  This function computes the Dirichlet-Voronoi tessellation of \code{X}
  (see \code{\link{dirichlet}})
  and extracts the vertices of the tiles of the tessellation.
  These vertices are the locations which are locally farthest away
  from \code{X}, that is, where the distance function
  of \code{X} reaches a local maximum.

  The Dirichlet tessellation of \code{X} is computed
  using \code{\link{dirichlet}}. The vertices of all
  tiles of the tessellation are extracted. Any vertex which
  lies on the boundary of the window of \code{X} is deleted.
  The remaining vertices are returned, as a point pattern,
  without duplicated entries.
}
\value{
  Point pattern (object of class \code{"ppp"}) in the same window
  as \code{X}.
}
\seealso{
 \code{\link{dirichlet}}
}
\examples{
  plot(dirichlet(cells))
  plot(dirichlet.vertices(cells), add=TRUE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{math}
\keyword{manip}
