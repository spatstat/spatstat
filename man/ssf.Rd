\name{ssf}
\alias{ssf}
\title{
  Spatially Sampled Function
}
\description{
  Create an object that represents a spatial function
  which has been evaluated or sampled at an irregular set of points.
}
\usage{
  ssf(loc, val)
}
\arguments{
  \item{loc}{
    The spatial locations at which the function has been evaluated.
    A point pattern (object of class \code{"ppp"}).
  }
  \item{val}{
    The function values at these locations.
    A numeric vector with one entry for each point of \code{loc},
    or a data frame with one row for each point of \code{loc}.
  }
}
\details{
  An object of class \code{"ssf"} represents a real-valued or
  vector-valued function that has been evaluated or sampled at an
  irregular set of points. An example would be a spatial covariate
  that has only been measured at certain locations.

  An object of this class also inherits the class \code{"ppp"},
  and is essentially the same as a marked point pattern, except
  for the class membership which enables it to be 
  handled in a different way.

  There are methods for \code{plot}, \code{print} etc; see
  \code{\link{plot.ssf}} and \code{\link{methods.ssf}}.

  Use \code{\link[spatstat]{unmark}} to extract only the point
  locations, and \code{\link{marks.ssf}} to extract only the function values.
}
\value{
  Object of class \code{"ssf"}.
}
\author{
  \adrian
}
\seealso{
  \code{\link{plot.ssf}}, 
  \code{\link{methods.ssf}}, 
  \code{\link{Smooth.ssf}}, 
  \code{\link{with.ssf}},
  \code{\link{[.ssf}}.
}
\examples{
  ssf(cells, nndist(cells, k=1:3))
}
\keyword{spatial}
\keyword{datagen}

