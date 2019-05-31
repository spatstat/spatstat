\name{angles.psp}
\alias{angles.psp}
\title{Orientation Angles of Line Segments}
\description{
  Computes the orientation angle of each line segment
  in a line segment pattern.
}
\usage{
  angles.psp(x, directed=FALSE)
}
\arguments{
  \item{x}{
    A line segment pattern (object of class \code{"psp"}).
  }
  \item{directed}{
    Logical flag. See details.
  }
}
\value{
  Numeric vector.
}
\details{
  For each line segment,  the angle of inclination to the \eqn{x}-axis
  (in radians) is computed,
  and the angles are returned as a numeric vector.

  If \code{directed=TRUE}, the directed angle of orientation
  is computed. The angle respects the
  sense of direction from \code{(x0,y0)} to \code{(x1,y1)}.
  The values returned are angles in the full range from \eqn{-\pi}{-\pi}
  to \eqn{\pi}{\pi}. The angle is computed as 
  \code{atan2(y1-y0,x1-x0)}. See \code{\link{atan2}}.

  If \code{directed=FALSE}, the undirected angle of orientation
  is computed. Angles differing by \eqn{\pi} are
  regarded as equivalent. The values returned are angles
  in the range from \eqn{0} to \eqn{\pi}{\pi}. These angles are
  computed by first computing the directed angle,
  then adding \eqn{\pi}{\pi} to any negative angles.
}
\seealso{
  \code{\link{psp}},
  \code{\link{marks.psp}},
  \code{\link{summary.psp}},
  \code{\link{midpoints.psp}},
  \code{\link{lengths.psp}},
  \code{\link{endpoints.psp}},
  \code{\link{extrapolate.psp}}.
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  b <- angles.psp(a)   
}
\author{
  \adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}

