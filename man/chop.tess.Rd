\name{chop.tess}
\alias{chop.tess}
\title{Subdivide a Window or Tessellation using a Set of Lines}
\description{
  Divide a given window into tiles
  delineated by a set of infinite straight lines, obtaining
  a tessellation of the window.
  Alternatively, given a tessellation, divide each tile of the
  tessellation into sub-tiles delineated by the lines.
}
\usage{
chop.tess(X, L)
}
\arguments{
  \item{X}{
    A window (object of class \code{"owin"}) or tessellation
    (object of class \code{"tess"}) to be subdivided by lines.
  }
  \item{L}{
    A set of infinite straight lines (object of class \code{"infline"})
  }
}
\details{
  The argument \code{L} should be a set of infinite straight lines in the plane
  (stored in an object \code{L} of class \code{"infline"} created by the
  function \code{\link{infline}}).
  
  If \code{X} is a window, then it is divided into tiles
  delineated by the lines in \code{L}.

  If \code{X} is a tessellation, then each tile of \code{X} is
  subdivided into sub-tiles delineated by the lines in \code{L}.

  The result is a tessellation. 
}
\section{Warning}{
  If \code{X} is a non-convex window, or a tessellation containing
  non-convex tiles, then \code{chop.tess(X,L)} may contain a tile
  which consists of several unconnected pieces. 
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{infline}},
  \code{\link{clip.infline}}
}
\examples{
  L <- infline(p=1:3, theta=pi/4)
  W <- square(4)
  chop.tess(W, L)
}
\keyword{spatial}
\keyword{math}
