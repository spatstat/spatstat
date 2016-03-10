\name{connected}  %DontDeclareMethods
\Rdversion{1.1}
\alias{connected}
\alias{connected.im}
\alias{connected.owin}
\title{
  Connected components
}
\description{
  Finds the topologically-connected components of a spatial object,
  such as the connected clumps of pixels in a binary image.
}
\usage{
connected(X, \dots)

\method{connected}{owin}(X, \dots, method="C")

\method{connected}{im}(X, \dots, background = NA, method="C")
}
\arguments{
  \item{X}{
    A spatial object such as a pixel image (object of class \code{"im"})
    or a window (object of class \code{"owin"}).
  }
  \item{background}{
    Optional. Treat pixels with this value 
    as being part of the background.
  }
  \item{method}{
    String indicating the algorithm to be used. Either \code{"C"}
    or \code{"interpreted"}. See Details.
  }
  \item{\dots}{
    Other arguments, not recognised by these methods.
  }
}
\details{
  The function \code{connected} is generic, with methods
  for pixel images (class \code{"im"}) and windows (class \code{"owin"})
  described here. There is also a method for point patterns
  described in \code{\link{connected.ppp}}.
  
  The functions described here compute the connected component transform
  (Rosenfeld and Pfalz, 1966)
  of a binary image or binary mask. The argument \code{X} is first
  converted into a pixel image with logical values. Then the algorithm
  identifies the connected components (topologically-connected clumps
  of pixels) in the foreground.

  Two pixels belong to the same connected component if they have the value
  \code{TRUE} and if they are neighbours (in the 8-connected
  sense). This rule is applied repeatedly until it terminates.
  Then each connected component
  contains all the pixels that can be reached by stepping from neighbour
  to neighbour.

  If \code{method="C"}, the computation is performed by a compiled C language
  implementation of the classical algorithm of Rosenfeld and Pfalz
  (1966). If \code{method="interpreted"}, the computation is performed
  by an \R implementation of the algorithm of Park et al (2000). 

  The result is a factor-valued image, with levels that correspond to
  the connected components. The Examples show how to extract each
  connected component as a separate window object.
}
\value{
  A pixel image (object of class \code{"im"}) with factor values.
  The levels of the factor correspond to the connected components.
}
\references{
  Park, J.-M., Looney, C.G. and Chen, H.-C. (2000)
  Fast connected component labeling algorithm using a divide and conquer
  technique. Pages 373-376 in
  S.Y. Shin (ed) \emph{Computers and Their Applications:} Proceedings of
  the ISCA 15th International Conference on Computers and Their
  Applications, March 29-31, 2000, New Orleans, Louisiana USA. ISCA
  2000, ISBN 1-880843-32-3. 

  Rosenfeld, A. and Pfalz, J.L. (1966)
  Sequential operations in digital processing.
  \emph{Journal of the Association for Computing Machinery} \bold{13}
  471-494.
}
\seealso{
  \code{\link{connected.ppp}},
  \code{\link{im.object}}, 
  \code{\link{tess}}
}
\section{Warnings}{
  It may be hard to distinguish different components 
  in the default plot because the colours of nearby components may be
  very similar. See the Examples for a randomised colour map.
  
  The algorithm for \code{method="interpreted"}
  can be very slow for large images (or images where
  the connected components include a large number of pixels).
}
\examples{
  data(cells)
  d <- distmap(cells, dimyx=256)
  X <- levelset(d, 0.06)
  plot(X)
  Z <- connected(X)
  plot(Z)

  # number of components
  nc <- length(levels(Z))
  # plot with randomised colour map
  plot(Z, col=hsv(h=sample(seq(0,1,length=nc), nc)))

  # how to extract the components as a list of windows
  W <- tiles(tess(image=Z))
}
\author{
  Original \R code by Julian Burgos, University of Washington.
  Adapted for \pkg{spatstat} by
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
