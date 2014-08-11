\name{opening} %DontDeclareMethods
\alias{opening}
\alias{opening.owin}
\alias{opening.psp}
\alias{opening.ppp}
\title{Morphological Opening}
\description{
  Perform morphological opening of a window, a line segment pattern
  or a point pattern.
}
\usage{
 opening(w, r, \dots)

 \method{opening}{owin}(w, r, \dots, polygonal=NULL)

 \method{opening}{ppp}(w, r, \dots)

 \method{opening}{psp}(w, r, \dots)
}
\arguments{
  \item{w}{
    A window (object of class \code{"owin"}
    or a line segment pattern (object of class \code{"psp"})
    or a point pattern (object of class \code{"ppp"}).
  }
  \item{r}{positive number: the radius of the opening.}
  \item{\dots}{
    extra arguments passed to \code{\link{as.mask}}
    controlling the pixel resolution, if a pixel approximation is used
  }
  \item{polygonal}{
    Logical flag indicating whether to compute a polygonal
    approximation to the erosion (\code{polygonal=TRUE}) or
    a pixel grid approximation (\code{polygonal=FALSE}).
  }
}
\value{
  If \code{r > 0}, an object of class \code{"owin"} representing the
  opened region. If \code{r=0}, the result is identical to \code{w}.
}
\details{
  The morphological opening (Serra, 1982)
  of a set \eqn{W} by a distance \eqn{r > 0}
  is the subset of points in \eqn{W} that can be 
  separated from the boundary of \eqn{W} by a circle of radius \eqn{r}.
  That is, a point \eqn{x} belongs to the opening 
  if it is possible to draw a circle of radius \eqn{r} (not necessarily
  centred on \eqn{x}) that has \eqn{x} on the inside
  and the boundary of \eqn{W} on the outside.
  The opened set is a subset of \code{W}.

  For a small radius \eqn{r}, the opening operation
  has the effect of smoothing out irregularities in the boundary of
  \eqn{W}. For larger radii, the opening operation removes promontories
  in the boundary. For very large radii, the opened set is empty.

  The algorithm applies \code{\link{erosion}} followed by
  \code{\link{dilation}}.
}
\seealso{
  \code{\link{closing}} for the opposite operation.

  \code{\link{dilation}}, \code{\link{erosion}} for the basic
  operations.  
  
  \code{\link{owin}},
  \code{\link{as.owin}} for information about windows.
}
\examples{
  v <- opening(letterR, 0.3, polygonal=interactive())
  plot(letterR, type="n", main="opening")
  plot(v, add=TRUE, col="grey")
  plot(letterR, add=TRUE)
}
\references{
  Serra, J. (1982)
  Image analysis and mathematical morphology.
  Academic Press.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
