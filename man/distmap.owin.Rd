\name{distmap.owin}
\alias{distmap.owin}
\title{Distance Map of Window}
\description{
  Computes the distance from each pixel to the nearest point
  in the given window.
}
\usage{
  \method{distmap}{owin}(X, \dots, discretise=FALSE, invert=FALSE)
}
\arguments{
  \item{X}{
    A window (object of class \code{"owin"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}}
    to control pixel resolution.
  }
  \item{discretise}{
    Logical flag controlling the choice of algorithm when \code{X} is
    a polygonal window. See Details.
  }
  \item{invert}{
    If \code{TRUE}, compute the distance transform of the
    complement of the window.
  }
}
\value{
  A pixel image (object of class \code{"im"}) whose greyscale values
  are the values of the distance map.
  The return value has an attribute \code{"bdry"}
  which is a pixel image.
}
\details{
  The ``distance map'' of a window \eqn{W} is the function
  \eqn{f} whose value \code{f(u)} is defined for any two-dimensional
  location \eqn{u} as the shortest distance from \eqn{u} to \eqn{W}.

  This function computes the distance map of the window \code{X}
  and returns the distance map as a pixel image. The greyscale value
  at a pixel \eqn{u} equals the distance from \eqn{u} to the nearest
  pixel in \code{X}.

  Additionally, the return value 
  has an attribute \code{"bdry"} which is 
  also a pixel image. The grey values in \code{"bdry"} give the
  distance from each pixel to the bounding rectangle of the image.

  If \code{X} is a binary pixel mask,
  the distance values computed are not the
  usual Euclidean distances. Instead the distance between two pixels
  is measured by the length of the
  shortest path connecting the two pixels. A path is a series of steps
  between neighbouring pixels (each pixel has 8 neighbours). 
  This is the standard `distance transform' algorithm of image
  processing (Rosenfeld and Kak, 1968; Borgefors, 1986).

  If \code{X} is a polygonal window, then exact Euclidean distances
  will be computed if \code{discretise=FALSE}. If \code{discretise=TRUE}
  then the window will first be converted to a binary pixel mask
  and the discrete path distances will be computed.

  The arguments \code{\dots} are passed to \code{\link{as.mask}}
  to control the pixel resolution.
  
  This function is a method for the generic \code{\link{distmap}}.
}
\seealso{
  \code{\link{distmap}},
  \code{\link{distmap.ppp}},
  \code{\link{distmap.psp}}
}
\examples{
  data(letterR)
  U <- distmap(letterR)
  \dontrun{
  plot(U)
  plot(attr(U, "bdry"))
  }
}
\references{
  Borgefors, G.
  Distance transformations in digital images.
  \emph{Computer Vision, Graphics and Image Processing} \bold{34}
  (1986) 344--371.
  
  Rosenfeld, A. and Pfalz, J.L.
  Distance functions on digital pictures.
  \emph{Pattern Recognition} \bold{1} (1968) 33-61.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
