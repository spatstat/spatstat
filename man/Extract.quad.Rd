\name{Extract.quad}
\alias{[.quad}
\title{Subset of Quadrature Scheme}
\description{
  Extract a subset of a quadrature scheme.
}
\usage{
  \method{[}{quad}(x, ...)
}
\arguments{
  \item{x}{
    A quadrature scheme (object of class \code{"quad"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{[.ppp}} to determine the
    subset.
  }
}
\value{
  A quadrature scheme (object of class \code{"quad"}).
}
\details{
  This function extracts a designated subset of a quadrature scheme.

  The function \code{[.quad} is a method for \code{\link{[}} for the
  class \code{"quad"}. It extracts a designated subset of a quadrature
  scheme.

  The subset to be extracted is determined by the arguments \code{\dots} which
  are interpreted by \code{\link{[.ppp}}. Thus it is possible to take the
  subset consisting of all quadrature points that lie inside a
  given region, or a subset of quadrature points identified by
  numeric indices.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{[.ppp}}.
}
\examples{
 data(nztrees)
 Q <- quadscheme(nztrees)
 W <- owin(c(0,148),c(0,95)) # a subwindow
 Q[W]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
