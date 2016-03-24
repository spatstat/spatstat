\name{erosion}  %DontDeclareMethods
\alias{erosion}
\alias{erosion.owin}
\alias{erosion.ppp}
\alias{erosion.psp}
\title{Morphological Erosion by a Disc}
\description{
  Perform morphological erosion of a window, a line segment pattern
  or a point pattern by a disc.
}
\usage{
 erosion(w, r, \dots)
 \method{erosion}{owin}(w, r, shrink.frame=TRUE, \dots,
   strict=FALSE, polygonal=NULL)
 \method{erosion}{ppp}(w, r,\dots)
 \method{erosion}{psp}(w, r,\dots)
}
\arguments{
  \item{w}{
    A window (object of class \code{"owin"}
    or a line segment pattern (object of class \code{"psp"})
    or a point pattern (object of class \code{"ppp"}).
  }
  \item{r}{positive number: the radius of erosion.}
  \item{shrink.frame}{logical: if \code{TRUE}, erode the bounding
    rectangle as well.}
  \item{\dots}{extra arguments to \code{\link{as.mask}}
    controlling the pixel resolution, if pixel approximation is used.}
  \item{strict}{Logical flag determining the fate of boundary pixels,
    if pixel approximation is used. See details.}
  \item{polygonal}{
    Logical flag indicating whether to compute a polygonal
    approximation to the erosion (\code{polygonal=TRUE}) or
    a pixel grid approximation (\code{polygonal=FALSE}).
  }
}
\value{
  If \code{r > 0}, an object of class \code{"owin"} representing the
  eroded region (or \code{NULL} if this region is empty).
  If \code{r=0}, the result is identical to \code{w}.
}
\details{
  The morphological erosion of a set \eqn{W} by a distance \eqn{r > 0}
  is the subset 
  consisting of all points \eqn{x \in W}{x in W} such that the
  distance from \eqn{x} to the boundary of \eqn{W} is greater than
  or equal to \eqn{r}. In other words it is the result of trimming
  a margin of width \eqn{r} off the set \eqn{W}.

  If \code{polygonal=TRUE} then a polygonal approximation
  to the erosion is computed.
  If \code{polygonal=FALSE} then a pixel approximation
  to the erosion is computed from the distance map of \code{w}.
  The arguments \code{"\dots"} are passed to \code{\link{as.mask}}
  to control the pixel resolution.
  The erosion consists of all pixels whose distance
  from the boundary of \code{w} is strictly greater than \code{r} (if
  \code{strict=TRUE}) or is greater than or equal to \code{r} (if
  \code{strict=FALSE}).
  
  When \code{w} is a window, the default (when \code{polygonal=NULL})
  is to compute a polygonal approximation if
  \code{w} is a rectangle or polygonal window, and to compute a
  pixel approximation if \code{w} is a window of type \code{"mask"}.

  If \code{shrink.frame} is false, the resulting window is given the
  same outer, bounding rectangle as the original window \code{w}.
  If \code{shrink.frame} is true, the original bounding rectangle
  is also eroded by the same distance \code{r}.

  To simply compute the area of the eroded window,
  use \code{\link{eroded.areas}}.
}
\seealso{
  \code{\link{dilation}} for the opposite operation.

  \code{\link{erosionAny}} for morphological erosion using any shape.
  
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{eroded.areas}}
}
\examples{
  plot(letterR, main="erosion(letterR, 0.2)")
  plot(erosion(letterR, 0.2), add=TRUE, col="red")
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{math}
