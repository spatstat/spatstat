\name{nobjects}
\alias{nobjects}
\alias{nobjects.ppp}
\alias{nobjects.ppx}
\alias{nobjects.psp}
\alias{nobjects.tess}
\title{
  Count Number of Geometrical Objects in a Spatial Dataset
}
\description{
  A generic function to count the number of geometrical objects in a
  spatial dataset.
}
\usage{
  nobjects(x)

  \method{nobjects}{ppp}(x)

  \method{nobjects}{ppx}(x)

  \method{nobjects}{psp}(x)

  \method{nobjects}{tess}(x)
}
\arguments{
  \item{x}{A dataset.}
}
\details{
  The generic function \code{nobjects} counts the number of
  geometrical objects in the spatial dataset \code{x}.

  The methods for point patterns (classes \code{"ppp"} and \code{"ppx"},
  embracing \code{"pp3"} and \code{"lpp"})
  count the number of points in the pattern.

  The method for line segment patterns (class \code{"psp"})
  counts the number of line segments in the pattern.

  The method for tessellations (class \code{"tess"})
  counts the number of tiles of the tessellation.
}
\value{
  A single integer.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{npoints}}
}
\examples{
   nobjects(redwood)
   nobjects(edges(letterR))
   nobjects(dirichlet(cells))
}
\keyword{spatial}
\keyword{manip}
