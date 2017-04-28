\name{as.linim} 
\alias{as.linim}
\alias{as.linim.linim}
\alias{as.linim.linfun}
\alias{as.linim.default}
\title{Convert to Pixel Image on Linear Network}
\description{
  Converts various kinds of data to a pixel image on a linear network.
}
\usage{
  as.linim(X, \dots)

  \method{as.linim}{linim}(X, \dots)

  \method{as.linim}{default}(X, L, \dots,
                             eps = NULL, dimyx = NULL, xy = NULL,
                             delta=NULL)

  \method{as.linim}{linfun}(X, L=domain(X), \dots,
                            eps = NULL, dimyx = NULL, xy = NULL,
                            delta=NULL)
}
\arguments{
  \item{X}{
    Data to be converted to a pixel image on a linear network.
  }
  \item{L}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{Additional arguments passed to \code{X} when \code{X}
    is a function.
  }
  \item{eps,dimyx,xy}{
    Optional arguments passed to \code{\link{as.mask}} to control
    the pixel resolution.
  }
  \item{delta}{
    Optional. Numeric value giving the approximate distance
    (in coordinate units) between successive
    sample points along each segment of the network.
  }
}
\details{
  This function converts the data \code{X} into a pixel image
  on a linear network, an object of class \code{"linim"}
  (see \code{\link{linim}}). 

  The argument \code{X} may be any of the following:
  \itemize{
    \item
    a function on a linear network, an object of class \code{"linfun"}.
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
  The conversion is performed by \code{\link{as.im}}.
  The arguments \code{eps}, \code{dimyx} and \code{xy}
  determine the pixel resolution.

  Next \code{Y} is converted to a pixel image on a linear network
  using \code{\link{linim}}. The argument \code{L} determines the
  linear network. If \code{L} is missing or \code{NULL},
  then \code{X} should be an object of class \code{"linim"},
  and \code{L} defaults to the linear network on which \code{X} is defined.

  In addition to converting the
  function to a pixel image, the algorithm also generates a fine grid of
  sample points evenly spaced along each segment of the network
  (with spacing at most \code{delta} coordinate units). The function values
  at these sample points are stored in the resulting object as a data frame
  (the argument \code{df} of \code{\link{linim}}). This mechanism allows
  greater accuracy for some calculations (such as
  \code{\link{integral.linim}}).
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
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
