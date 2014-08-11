\name{Extract.influence.ppm}
\alias{[.influence.ppm}
\title{Extract Subset of Influence Object}
\description{
  Extract a subset of an influence object, 
  or extract the influence values at specified locations.
}
\usage{
  \method{[}{influence.ppm}(x, i, ...)
}
\arguments{
  \item{x}{
    A influence object (of class \code{"influence.ppm"})
    computed by \code{\link{influence.ppm}}.
  }
  \item{i}{
    Subset index (passed to \code{\link{[.ppp}}).
    Either a spatial window (object of class \code{"owin"})
    or an integer index.
  }
  \item{\dots}{
     Ignored.
  }
}
\value{
  Another object of class \code{"influence.ppm"}.
}
\details{
  An object of class \code{"influence.ppm"} contains the values of the
  likelihood influence for a point process model, computed by
  \code{\link{influence.ppm}}. This is effectively a marked point
  pattern obtained by marking each of the original data points
  with its likelihood influence.
  
  This function extracts a designated subset of the influence values,
  either as another influence object, or as a vector of numeric values.

  The function \code{[.influence.ppm} is a method for \code{\link{[}} for the
  class \code{"influence.ppm"}. The argument \code{i} should be
  an index applicable to a point pattern. It may be either
  a spatial window (object of class \code{"owin"})
  or a sequence index.
  The result will be another influence object
  (of class \code{influence.ppm}).

  To extract the influence values as a numeric vector,
  use \code{marks(as.ppp(x))}.
}
\seealso{
  \code{\link{influence.ppm}}.
}
\examples{
   fit <- ppm(cells, ~x)
   infl <- influence(fit)
   b <- owin(c(0.1, 0.3), c(0.2, 0.4))
   infl[b]
   infl[1:5]
   marks(as.ppp(infl))[1:3]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
