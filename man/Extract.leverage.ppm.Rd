\name{Extract.leverage.ppm}
\alias{[.leverage.ppm}
\title{Extract Subset of Leverage Object}
\description{
  Extract a subset of a leverage map, 
  or extract the leverage values at specified locations.
}
\usage{
  \method{[}{leverage.ppm}(x, i, ...)
}
\arguments{
  \item{x}{
    A leverage object (of class \code{"leverage.ppm"})
    computed by \code{\link{leverage.ppm}}.
  }
  \item{i}{
    Subset index (passed to \code{\link{[.im}}).
    Either a spatial window (object of class \code{"owin"})
    or a spatial point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Further arguments passed to \code{\link{[.im}},
    especially the argument \code{drop}.
  }
}
\value{
  Another object of class \code{"leverage.ppm"},
  or a vector of numeric values of leverage.
}
\details{
  An object of class \code{"leverage.ppm"} contains the values of the
  leverage function for a point process model, computed by
  \code{\link{leverage.ppm}}.

  This function extracts a designated subset of the leverage values,
  either as another leverage object, or as a vector of numeric values.

  The function \code{[.leverage.ppm} is a method for \code{\link{[}} for the
  class \code{"leverage.ppm"}. The argument \code{i} should be either
  \itemize{
    \item 
     a spatial window (object of class \code{"owin"})
     determining a region where the leverage map is required.
     The result will typically be another leverage map
     (object of class \code{leverage.ppm}).
   \item 
     a spatial point pattern (object of class \code{"ppp"})
     specifying locations at which the leverage values are required.
     The result will be a numeric vector.
  }
  The subset operator for images, \code{\link{[.im}}, is applied to
  the leverage map. If this yields a pixel image, then the result of
  \code{\link{[.leverage.ppm}} is another leverage object. Otherwise,
  a vector containing the numeric values of leverage is returned.
}
\seealso{
  \code{\link{leverage.ppm}}.
}
\examples{
   fit <- ppm(cells, ~x)
   lev <- leverage(fit)
   b <- owin(c(0.1, 0.3), c(0.2, 0.4))
   lev[b]
   lev[cells]
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{manip}
