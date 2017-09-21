\name{as.function.leverage.ppm}
\alias{as.function.leverage.ppm}
\title{
  Convert Leverage Object to Function of Coordinates
}
\description{
  Converts an object of class \code{"leverage.ppm"} to a function of the
  \eqn{x} and \eqn{y} coordinates.
}
\usage{
 \method{as.function}{leverage.ppm}(x, ...)
}
\arguments{
  \item{x}{
    Object of class \code{"leverage.ppm"}
    produced by \code{\link{leverage.ppm}}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  An object of class \code{"leverage.ppm"} represents the leverage
  function of a fitted point process model. This command converts the object
  to a \code{function(x,y)}
  where the arguments \code{x} and \code{y} are (vectors of) spatial
  coordinates. This function returns the leverage values at the
  specified locations (calculated by referring to the nearest location
  where the leverage has been computed).
}
\value{
  A function in the \R language, also belonging to the
  class \code{"funxy"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{as.im.leverage.ppm}}
}
\examples{
  X <- rpoispp(function(x,y) { exp(3+3*x) })
  fit <- ppm(X ~x+y)
  lev <- leverage(fit)
  f <- as.function(lev)
  
  f(0.2, 0.3)  # evaluate at (x,y) coordinates
  y <- f(X)    # evaluate at a point pattern
}
\keyword{spatial}
\keyword{manip}
