\name{as.linim}
\alias{as.linim}
\alias{as.linim.linim}
\alias{as.linim.default}
\title{Convert to Pixel Image on Linear Network}
\description{
  Converts various kinds of data to a pixel image on a linear network.
}
\usage{
  as.linim(X, \dots)

  \method{as.linim}{linim}(X, \dots)

  \method{as.linim}{default}(X, L, \dots)
}
\arguments{
  \item{X}{
    Data to be converted to a pixel image on a linear network.
  }
  \item{L}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{Additional arguments passed to \code{X} when \code{X}
    is a function, or arguments passed to \code{\link{as.mask}}
    which determine the pixel array geometry.
  }
}
\details{
  This function converts the data \code{X} into a pixel image
  on a linear network, an object of class \code{"linim"}
  (see \code{\link{linim}}). 

  The argument \code{X} may be any of the following:
  \itemize{
    \item
    a pixel image on a linear network, an object of class
    \code{"linim"}.
    \item
    a pixel image, an object of class \code{"im"}.
    \item
    any type of data acceptable to \code{\link{as.im}},
    such as a function, numeric value, or window.
  }

  First \code{X} is converted to a pixel image object \code{Y}
  (object of class \code{"im"}).
  The conversion is performed by \code{\link{as.im}}
  using any relevant arguments specified in \code{\dots}. For example
  the argument \code{eps} could be used to change the pixel resolution.

  Next \code{Y} is converted to a pixel image on a linear network
  using \code{\link{linim}}. The argument \code{L} determines the
  linear network. If \code{L} is missing or \code{NULL},
  then \code{X} should be an object of class \code{"linim"},
  and \code{L} defaults to the linear network on which \code{X} is defined.
}
\value{
  An image object on a linear network; an object of class \code{"linim"}.
}
\seealso{
  \code{\link{as.im}}
}
\examples{
  f <- function(x,y){ x + y }
  plot(as.linim(f, simplenet))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
