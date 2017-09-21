\name{as.function.owin}
\alias{as.function.owin}
\title{
  Convert Window to Indicator Function
}
\description{
  Converts a spatial window to a function of the
  \eqn{x} and \eqn{y} coordinates returning the value 1 inside the
  window and 0 outside.
}
\usage{
 \method{as.function}{owin}(x, \dots)
}
\arguments{
  \item{x}{
    Pixel image (object of class \code{"owin"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This command converts a spatial window
  (object of class \code{"owin"}) to a \code{function(x,y)}
  where the arguments \code{x} and \code{y} are (vectors of) spatial
  coordinates. This is the indicator function of the window:
  it returns the value 1 for locations inside the window,
  and returns 0 for values outside the window.
}
\value{
  A function in the \R language.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{as.im.owin}}
}
\examples{
  W <- Window(humberside)
  f <- as.function(W)
  f(5000, 4500)
  f(123456, 78910)
  X <- runifpoint(5, Frame(humberside))
  f(X)
}
\keyword{spatial}
\keyword{manip}
