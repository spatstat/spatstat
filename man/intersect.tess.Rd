\name{intersect.tess}
\alias{intersect.tess}
\title{Intersection of Two Tessellations}
\description{
  Yields the intersection of two tessellations,
  or the intersection of a tessellation with a window.
}
\usage{
   intersect.tess(X, Y, \dots, keepmarks=FALSE, sep="x")
}
\arguments{
  \item{X,Y}{Two tessellations (objects of class \code{"tess"}),
    or windows (objects of class \code{"tess"}),
    or other data that can be converted to 
    tessellations by \code{\link{as.tess}}.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}}
    to control the discretisation, if required.
  }
  \item{keepmarks}{
    Logical value. If \code{TRUE}, the marks attached to the
    tiles of \code{X} and \code{Y} will be retained as marks of
    the intersection tiles. 
  }
  \item{sep}{
    Character string used to separate the names of tiles from \code{X}
    and from \code{Y}, when forming the name of the tiles of the
    intersection.
  }
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\details{
  A tessellation is a collection of disjoint spatial regions
  (called \emph{tiles}) that fit together to form a larger spatial
  region. See \code{\link{tess}}.

  If \code{X} and \code{Y} are not tessellations, they are first
  converted into tessellations by \code{\link{as.tess}}.

  The function \code{intersect.tess} then computes the intersection between
  the two tessellations. This is another tessellation, each of whose
  tiles is the intersection of a tile from \code{X} and a tile from \code{Y}.

  One possible use of this function is to slice a window \code{W} into
  subwindows determined by a tessellation. See the Examples.
}

\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{tess}},
  \code{\link{as.tess}},
  \code{\link{intersect.owin}}
}
\examples{
  opa <- par(mfrow=c(1,3))
# polygon
  data(letterR)
  plot(letterR)
# tessellation of rectangles
  X <- tess(xgrid=seq(2, 4, length=10), ygrid=seq(0, 3.5, length=8))
  plot(X)
  plot(intersect.tess(X, letterR))

  A <- runifpoint(10)
  B <- runifpoint(10)
  plot(DA <- dirichlet(A))
  plot(DB <- dirichlet(B))
  plot(intersect.tess(DA, DB))
  par(opa)

  marks(DA) <- 1:10
  marks(DB) <- 1:10
  plot(Z <- intersect.tess(DA,DB, keepmarks=TRUE))
  mZ <- marks(Z)
  tZ <- tiles(Z)
  for(i in which(mZ[,1] == 3)) plot(tZ[[i]], add=TRUE, col="pink")
}
\keyword{spatial}
\keyword{math}
