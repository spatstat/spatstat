\name{dirichlet.vertices}
\alias{dirichlet.vertices}
\alias{dirichlet.edges}
\title{
   Vertices and Edges of Dirichlet Tessellation
}
\description{
  Computes the Dirichlet-Voronoi tessellation of a point pattern
  and extracts the vertices or edges of the tiles.
}
\usage{
dirichlet.vertices(X)

dirichlet.edges(X)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
}
\details{
  These function compute the Dirichlet-Voronoi tessellation of \code{X}
  (see \code{\link{dirichlet}})
  and extract the vertices or edges of the tiles of the tessellation.

  The Dirichlet vertices are the spatial locations which are locally
  farthest away from \code{X}, that is, where the distance function
  of \code{X} reaches a local maximum.

  The Dirichlet edges are the dividing lines equally distant between
  a pair of points of \code{X}.

  The Dirichlet tessellation of \code{X} is computed
  using \code{\link{dirichlet}}. The vertices or edges of all
  tiles of the tessellation are extracted.

  For \code{dirichlet.vertices}, any vertex which
  lies on the boundary of the window of \code{X} is deleted.
  The remaining vertices are returned, as a point pattern,
  without duplicated entries.
}
\value{
  \code{dirichlet.vertices}
  returns a point pattern (object of class \code{"ppp"}) in the same window
  as \code{X}.

  \code{dirichlet.edges} returns a line segment pattern (object of
  class \code{"psp"}).
}
\seealso{
 \code{\link{dirichlet}}
}
\examples{
  plot(dirichlet(cells))

  plot(dirichlet.vertices(cells), add=TRUE)

  ed <- dirichlet.edges(cells)
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
