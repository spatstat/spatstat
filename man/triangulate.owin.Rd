\name{triangulate.owin}
\alias{triangulate.owin}
\title{
  Decompose Window into Triangles
}
\description{
  Given a spatial window, this function decomposes the window
  into disjoint triangles. 
  The result is a tessellation of the window
  in which each tile is a triangle.
}
\usage{
triangulate.owin(W)
}
\arguments{
  \item{W}{Window (object of class \code{"owin"}).}
}
\details{
  The window \code{W} will be decomposed into disjoint triangles.
  The result is a tessellation of \code{W} in which each tile is a
  triangle. All triangle vertices lie on the boundary
  of the original polygon.
  
  The window is first converted to a polygonal window using
  \code{\link{as.polygonal}}. The vertices of the polygonal window
  are extracted, and the Delaunay triangulation of these vertices
  is computed using \code{\link{delaunay}}. Each Delaunay triangle
  is intersected with the window: if the result is not a triangle, 
  the triangulation procedure is applied recursively to this smaller polygon. 
}
\value{
  Tessellation (object of class \code{"tess"}).
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{tess}},
  \code{\link{delaunay}},
  \code{\link{as.polygonal}}
}
\examples{
  plot(triangulate.owin(letterR))
}
\keyword{spatial}
\keyword{manip}

