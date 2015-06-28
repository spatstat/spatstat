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
  triangle.
  
  The window is first converted to a polygonal window using
  \code{\link{as.polygonal}}. The vertices of the polygonal window
  are extracted, and the Delaunay triangulation of these vertices
  is computed using \code{\link{delaunay}}. The algorithm then
  efficiently determines which triangles lie inside the window,
  and assembles these triangles into a tessellation.

  The vertices of the triangles are vertices of the original polygon.
}
\value{
  Tessellation (object of class \code{"tess"}).
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
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

