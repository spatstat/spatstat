\name{extrapolate.psp}
\alias{extrapolate.psp}
\title{
  Extrapolate Line Segments to Obtain Infinite Lines
}
\description{
  Given a spatial pattern of line segments,
  extrapolate the segments to infinite lines.
}
\usage{
extrapolate.psp(x, \dots)
}
\arguments{
  \item{x}{
    Spatial pattern of line segments (object of class \code{"psp"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  Each line segment in the pattern \code{x} is
  extrapolated to an infinite line, drawn through its two endpoints.
  The resulting pattern of infinite lines
  is returned as an object of class \code{"infline"}.

  If a segment's endpoints are identical
  (so that it has zero length)
  the resulting infinite line is vertical (i.e. parallel to the \eqn{y}
  coordinate axis).
}
\value{
  An object of class \code{"infline"} representing the pattern of
  infinite lines. See \code{\link{infline}} for details of
  structure.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{psp}},
  \code{\link{infline}}

  \code{\link{midpoints.psp}},
  \code{\link{lengths_psp}}
  \code{\link{angles.psp}},
  \code{\link{endpoints.psp}}.
}
\examples{
   X <- psp(runif(4), runif(4), runif(4), runif(4), window=owin())
   Y <- extrapolate.psp(X)
   plot(X, col=3, lwd=4)
   plot(Y, lty=3)
   Y
}
\keyword{spatial}
\keyword{math}
