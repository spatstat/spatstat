\name{dirichletVertices}
\alias{dirichletVertices}
\alias{dirichletEdges}
\title{
   Vertices and Edges of Dirichlet Tessellation
}
\description{
  Computes the Dirichlet-Voronoi tessellation of a point pattern
  and extracts the vertices or edges of the tiles.
}
\usage{
dirichletVertices(X)

dirichletEdges(X, clip=TRUE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{clip}{
    Logical value specifying whether to clip the tile edges to the
    window. See Details.
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

  For \code{dirichletVertices}, any vertex which
  lies on the boundary of the window of \code{X} is deleted.
  The remaining vertices are returned, as a point pattern,
  without duplicated entries.

  For \code{dirichletEdges}, the edges are initially computed
  inside the rectangle \code{Frame(X)}.
  Then if \code{clip=TRUE} (the default), these edges are intersected
  with \code{Window(X)}, which may cause an edge to be broken into
  several pieces.
}
\value{
  \code{dirichletVertices}
  returns a point pattern (object of class \code{"ppp"}) in the same window
  as \code{X}.

  \code{dirichletEdges} returns a line segment pattern (object of
  class \code{"psp"}).
}
\seealso{
 \code{\link{dirichlet}},
 \code{\link{dirichletAreas}}
}
\examples{
  plot(dirichlet(cells))

  plot(dirichletVertices(cells), add=TRUE)

  ed <- dirichletEdges(cells)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
\keyword{manip}
