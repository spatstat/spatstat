\name{quantess}
\alias{quantess}
\alias{quantess.owin}
\alias{quantess.ppp}
\alias{quantess.im}
\title{Quantile Tessellation}
\description{
  Divide space into tiles which contain equal amounts of stuff.
}
\usage{
quantess(M, Z, n, \dots)

\method{quantess}{owin}(M, Z, n, \dots, type=2, origin=c(0,0), eps=NULL)

\method{quantess}{ppp}(M, Z, n, \dots, type=2, origin=c(0,0), eps=NULL)

\method{quantess}{im}(M, Z, n, \dots, type=2, origin=c(0,0))
}
\arguments{
  \item{M}{
    A spatial object (such as a window, point pattern
    or pixel image) determining the weight or amount of stuff
    at each location.
  }
  \item{Z}{
    A spatial covariate (a pixel image or a \code{function(x,y)})
    or one of the strings \code{"x"} or \code{"y"} indicating the
    Cartesian coordinates \eqn{x} or \eqn{y}, or one of the strings
    \code{"rad"} or \code{"ang"} indicating polar coordinates.
    The range of values of \code{Z} will be broken into \code{n}
    bands containing equal amounts of stuff.
  }
  \item{n}{
    Number of bands. A positive integer.
  }
  \item{type}{
    Integer specifying the rule for calculating quantiles.
    Passed to \code{\link[stats]{quantile.default}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{quadrats}} or
    \code{\link{tess}} defining another tessellation
    which should be intersected with the quantile tessellation.
  }
  \item{origin}{
    Location of the origin of polar coordinates,
    if \code{Z="rad"} or \code{Z="ang"}.
    Either a numeric vector of length 2 giving the location,
    or a point pattern containing only one point,
    or a list with two entries named \code{x} and \code{y},
    or one of the character strings
    \code{"centroid"}, \code{"midpoint"},
    \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"},
    \code{"topleft"}, \code{"bottomleft"}, \code{"topright"} or
    \code{"bottomright"} (partially matched).
  }
  \item{eps}{
    Optional. The size of pixels in the approximation which is used
    to compute the quantiles. A positive numeric value, or vector of two
    positive numeric values.
  }
}
\details{
  A \emph{quantile tessellation} is a division of space into
  pieces which contain equal amounts of stuff.

  The function \code{quantess} 
  computes a quantile tessellation and 
  returns the tessellation itself.  
  The function \code{quantess} is generic, with methods for
  windows (class \code{"owin"}), point patterns (\code{"ppp"})
  and pixel images (\code{"im"}). 

  The first argument \code{M} (for mass) specifies the spatial
  distribution of stuff that is to be divided. If \code{M} is a window,
  the \emph{area} of the window is to be divided into \code{n} equal pieces.
  If \code{M} is a point pattern, the \emph{number of points} in the
  pattern is to be divided into \code{n} equal parts, as far as
  possible. If \code{M} is a pixel image, the pixel values are
  interpreted as weights, and the \emph{total weight} is to be divided
  into \code{n} equal parts.  

  The second argument
  \code{Z} is a spatial covariate. The range of values of \code{Z}
  will be divided into \code{n} bands, each containing
  the same total weight. That is, we determine the quantiles of \code{Z}
  with weights given by \code{M}.

  For convenience, additional arguments \code{\dots} can be given,
  to further subdivide the tiles of the tessellation. These arguments
  should be recognised by one of the functions
  \code{\link{quadrats}} or \code{\link{tess}}. The tessellation
  determined by these arguments is intersected with the quantile tessellation.

  The result of \code{quantess} is a tessellation of \code{as.owin(M)}
  determined by the quantiles of \code{Z}.
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\author{
  Original idea by Ute Hahn.
  Implemented in \code{spatstat} by
  \spatstatAuthors.
}
\seealso{
  \code{\link{tess}}, \code{\link{quadrats}}, \code{\link{quantile}},
  \code{\link{tilenames}}
}
\examples{
  plot(quantess(letterR, "x", 5))

  plot(quantess(bronzefilter, "x", 6))
  points(unmark(bronzefilter))

  plot(quantess(letterR, "rad", 7, origin=c(2.8, 1.5)))
  plot(quantess(letterR, "ang", 7, origin=c(2.8, 1.5)))

  opa <- par(mar=c(0,0,2,5))
  A <- quantess(Window(bei), bei.extra$elev, 4)
  plot(A, ribargs=list(las=1))
  
  B <- quantess(bei, bei.extra$elev, 4)
  tilenames(B) <- paste(spatstat.utils::ordinal(1:4), "quartile")
  plot(B, ribargs=list(las=1))
  points(bei, pch=".", cex=2, col="white")
  par(opa)
}
\keyword{spatial}
\keyword{manip}
