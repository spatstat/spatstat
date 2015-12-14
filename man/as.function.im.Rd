\name{as.function.im}
\alias{as.function.im}
\title{
  Convert Pixel Image to Function of Coordinates
}
\description{
  Converts a pixel image to a function of the
  \eqn{x} and \eqn{y} coordinates.
}
\usage{
 \method{as.function}{im}(x, ...)
}
\arguments{
  \item{x}{
    Pixel image (object of class \code{"im"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This command converts a pixel image
  (object of class \code{"im"}) to a \code{function(x,y)}
  where the arguments \code{x} and \code{y} are (vectors of) spatial
  coordinates. This function returns the pixel values at the
  specified locations.
}
\value{
  A function in the \R language, also belonging to the
  class \code{"funxy"}.
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
  \code{\link{[.im}}
}
\examples{
  d <- density(cells)
  f <- as.function(d)
  f(0.1, 0.3)
}
\keyword{spatial}
\keyword{manip}
