\name{is.convex}
\alias{is.convex}
\title{Test Whether a Window is Convex}
\description{
  Determines whether a window is convex.
}
\usage{
  is.convex(x)
}
\arguments{
  \item{x}{
    Window (object of class \code{"owin"}).
  }
}
\value{
  Logical value, equal to \code{TRUE} if \code{x} is convex.
}
\details{
  If \code{x} is a rectangle, the result is TRUE.

  If \code{x} is polygonal, the result is TRUE if \code{x} consists of a
  single polygon and this polygon is equal to the minimal convex hull
  of its vertices computed by \code{\link[grDevices]{chull}}. 

  If \code{x} is a mask, the algorithm first extracts all boundary
  pixels of \code{x} using \code{\link{vertices}}. Then it computes
  the (polygonal) convex hull \eqn{K} of the boundary pixels.
  The result is TRUE if every boundary pixel lies within 
  one pixel diameter of an edge of \eqn{K}.
}
\seealso{
  \code{\link{owin}},
  \code{\link{convexhull.xy}},
  \code{\link{vertices}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
