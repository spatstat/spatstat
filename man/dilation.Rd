\name{dilation} %DontDeclareMethods
\alias{dilation}
\alias{dilation.owin}
\alias{dilation.ppp}
\alias{dilation.psp}
\title{Morphological Dilation}
\description{
  Perform morphological dilation of a window, a line segment pattern
  or a point pattern
}
\usage{
 dilation(w, r, \dots)
 \method{dilation}{owin}(w, r, \dots, polygonal=NULL, tight=TRUE)
 \method{dilation}{ppp}(w, r, \dots, polygonal=TRUE, tight=TRUE)
 \method{dilation}{psp}(w, r, \dots, polygonal=TRUE, tight=TRUE)
}
\arguments{
  \item{w}{
    A window (object of class \code{"owin"}
    or a line segment pattern (object of class \code{"psp"})
    or a point pattern (object of class \code{"ppp"}).
  }
  \item{r}{positive number: the radius of dilation.}
  \item{\dots}{extra arguments passed to \code{\link{as.mask}}
    controlling the pixel resolution, if the pixel approximation is
    used.
  }
  \item{polygonal}{
    Logical flag indicating whether to compute a polygonal
    approximation to the dilation (\code{polygonal=TRUE}) or
    a pixel grid approximation (\code{polygonal=FALSE}).
  }
  \item{tight}{
    Logical flag indicating whether the bounding frame of the window
    should be taken as the smallest rectangle enclosing the dilated region
    (\code{tight=TRUE}), or should be the
    dilation of the bounding frame of \code{w} (\code{tight=FALSE}).
  }
}
\value{
  If \code{r > 0}, an object of class \code{"owin"} representing the
  dilated region. If \code{r=0}, the result is identical to \code{w}.
}
\details{
  The morphological dilation of a set \eqn{W} by a distance \eqn{r > 0}
  is the set consisting of all points lying at most \eqn{r} units
  away from \eqn{W}. Effectively, dilation adds a margin of width
  \eqn{r} onto the set \eqn{W}.

  If \code{polygonal=TRUE} then a polygonal approximation
  to the dilation is computed.
  If \code{polygonal=FALSE} then a pixel approximation
  to the dilation is computed from the distance map of \code{w}.
  The arguments \code{"\dots"} are passed to \code{\link{as.mask}}
  to control the pixel resolution.

  When \code{w} is a window, the default (when \code{polygonal=NULL})
  is to compute a polygonal approximation if
  \code{w} is a rectangle or polygonal window, and to compute a
  pixel approximation if \code{w} is a window of type \code{"mask"}.
}
\seealso{
  \code{\link{erosion}} for the opposite operation.
  
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  plot(dilation(letterR, 0.2))
  plot(letterR, add=TRUE, lwd=2, border="red")

  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  plot(dilation(X, 0.1))
  plot(X, add=TRUE, col="red")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
