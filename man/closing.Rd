\name{closing} %DontDeclareMethods
\alias{closing}
\alias{closing.owin}
\alias{closing.ppp}
\alias{closing.psp}
\title{Morphological Closing}
\description{
  Perform morphological closing of a window, a line segment pattern
  or a point pattern.
}
\usage{
 closing(w, r, \dots)

 \method{closing}{owin}(w, r, \dots, polygonal=NULL)

 \method{closing}{ppp}(w, r, \dots, polygonal=TRUE)

 \method{closing}{psp}(w, r, \dots, polygonal=TRUE)
}
\arguments{
  \item{w}{
    A window (object of class \code{"owin"}
    or a line segment pattern (object of class \code{"psp"})
    or a point pattern (object of class \code{"ppp"}).
  }
  \item{r}{positive number: the radius of the closing.}
  \item{\dots}{extra arguments passed to \code{\link{as.mask}}
    controlling the pixel resolution, if a pixel approximation is used}
  \item{polygonal}{
    Logical flag indicating whether to compute a polygonal
    approximation to the erosion (\code{polygonal=TRUE}) or
    a pixel grid approximation (\code{polygonal=FALSE}).
    Ignored if \code{gpclib} is disabled.
  }
}
\value{
  If \code{r > 0}, an object of class \code{"owin"} representing the
  closed region. If \code{r=0}, the result is identical to \code{w}.
}
\details{
  The morphological closing (Serra, 1982)
  of a set \eqn{W} by a distance \eqn{r > 0}
  is the set of all points that cannot be
  separated from \eqn{W} by any circle of radius \eqn{r}.
  That is, a point \eqn{x} belongs to the closing \eqn{W*}
  if it is impossible to draw any circle of radius \eqn{r} that
  has \eqn{x} on the inside and \eqn{W} on the outside.
  The closing \eqn{W*} contains the original set \eqn{W}.

  For a small radius \eqn{r}, the closing operation
  has the effect of smoothing out irregularities in the boundary of
  \eqn{W}. For larger radii, the closing operation smooths out
  concave features in the boundary. For very large radii,
  the closed set \eqn{W*} becomes more and more convex.

  The algorithm applies \code{\link{dilation}} followed by
  \code{\link{erosion}}. 

  Polygonal calculations require the \code{gpclib}
  package which is subject to licence restrictions.
  It is enabled by \code{spatstat.options(gpclib=TRUE)}.
  See \code{\link{licence.polygons}}.
}
\seealso{
  \code{\link{opening}} for the opposite operation.

  \code{\link{dilation}}, \code{\link{erosion}} for the basic
  operations.  
  
  \code{\link{owin}},
  \code{\link{as.owin}} for information about windows.
}
\examples{
  if(spatstat.options("gpclib")) {
    v <- closing(letterR, 0.25)
  } else {
    v <- opening(letterR, 0.25, dimyx=256, polygonal=FALSE)
  }
  plot(v, main="closing")
  plot(letterR, add=TRUE)
}
\references{
  Serra, J. (1982)
  Image analysis and mathematical morphology.
  Academic Press.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
