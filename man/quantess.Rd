\name{quantess}
\alias{quantsplit}
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

quantsplit(M, Z, n, \dots)

\method{quantess}{owin}(M, Z, n, \dots)

\method{quantess}{ppp}(M, Z, n, \dots)

\method{quantess}{im}(M, Z, n, \dots)
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
    \eqn{x} or \eqn{y} coordinate.
    The range of values of \code{Z} will be broken into \code{n}
    bands containing equal amounts of stuff.
  }
  \item{n}{
    Number of bands. A positive integer.
  }
  \item{\dots}{
    Ignored.
    Additional arguments passed to \code{\link{quadrats}} or
    \code{\link{tess}} defining another tessellation
    which should be intersected with the quantile tessellation.
  }
}
\details{
  A \emph{quantile tessellation} is a division of space into
  pieces which contain equal amounts of stuff.

  The functions \code{quantess} and \code{quantsplit}
  both compute a quantile tessellation;
  \code{quantess} returns the tessellation itself, while
  \code{quantsplit} splits the data \code{M} according to this tessellation.
  
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
  to further subdivide the tiles of the tessellation.

  The result of \code{quantess} is a tessellation of \code{as.owin(M)}
  determined by the quantiles of \code{Z}.
  The result of \code{quantsplit} is a list of spatial objects
  obtained by splitting \code{M} according to the quantile tessellation.
}
\value{
  The result of \code{quantess} is a
  tessellation (object of class \code{"tess"}).

  The result of \code{quantsplit} is a list (of class \code{"solist"})
  of objects of the same type as \code{M}.
}
\author{
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}.
Based on work by Ute Hahn.
}
\seealso{
  \code{\link{tess}}, \code{\link{quadrats}}, \code{\link{quantile}},
  \code{\link{tilenames}}
}
\examples{
  plot(quantess(letterR, "x", 5))

  plot(quantess(bronzefilter, "x", 6))
  points(unmark(bronzefilter))

  opa <- par(mar=c(0,0,2,5))
  A <- quantess(Window(bei), bei.extra$elev, 4)
  plot(A, ribargs=list(las=1))
  
  B <- quantess(bei, bei.extra$elev, 4)
  tilenames(B) <- paste(ordinal(1:4), "quartile")
  plot(B, ribargs=list(las=1))
  points(bei, pch=".", cex=2, col="white")
  par(opa)
}
\keyword{spatial}
\keyword{manip}
