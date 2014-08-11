\name{npoints}
\alias{npoints}
\alias{npoints.ppp}
\alias{npoints.pp3}
\alias{npoints.ppx}
\title{Number of Points in a Point Pattern}
\description{
  Returns the number of points in a point pattern of any kind.
}
\usage{
  npoints(x)
  \method{npoints}{ppp}(x)
  \method{npoints}{pp3}(x)
  \method{npoints}{ppx}(x)
}
\arguments{
  \item{x}{
    A point pattern (object of class \code{"ppp"},
    \code{"pp3"}, \code{"ppx"} or some other suitable class).
  }
}
\value{
  Integer.
}
\details{
  This function returns the number of points in
  a point pattern. The function \code{npoints} is generic
  with methods for the classes \code{"ppp"}, \code{"pp3"},
  \code{"ppx"} and possibly other classes.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{print.pp3}},
  \code{\link{print.ppx}}.
}
\examples{
   data(cells)
   npoints(cells)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
