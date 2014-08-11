\name{pairdist.ppp}
\alias{pairdist.ppp}
\title{Pairwise distances}
\description{
  Computes the matrix of distances between all pairs of points
  in a point pattern.
}
\usage{
  \method{pairdist}{ppp}(X, \dots, periodic=FALSE, method="C", squared=FALSE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{periodic}{
    Logical. Specifies whether to apply a periodic edge correction.
  }
  \item{method}{
    String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
    Usually not specified.
  }
  \item{squared}{
    Logical. If \code{squared=TRUE}, the squared distances are
    returned instead (this computation is faster).
  }
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance
  between the points numbered \code{i} and \code{j}.
}
\details{
  This is a method for the generic function \code{pairdist}.

  Given a point pattern \code{X} (an object of class \code{"ppp"}),
  this function computes the Euclidean distances between all pairs of
  points in \code{X}, and returns the matrix of distances.

  Alternatively if \code{periodic=TRUE} and the window containing \code{X} is a
  rectangle, then the distances will be computed in the `periodic'
  sense (also known as `torus' distance): opposite edges of the
  rectangle are regarded as equivalent.
  This is meaningless if the window is not a rectangle.

  If \code{squared=TRUE} then the \emph{squared} Euclidean distances
  \eqn{d^2} are returned, instead of the Euclidean distances \eqn{d}.
  The squared distances are faster to calculate, and are sufficient for
  many purposes (such as finding the nearest neighbour of a point).
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. The C code is somewhat faster.
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{pairdist.default}},
  \code{\link{pairdist.psp}},
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link{Kest}}
}
\examples{
   data(cells)
   d <- pairdist(cells)
   d <- pairdist(cells, periodic=TRUE)
   d <- pairdist(cells, squared=TRUE)
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and 
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}
\keyword{math}

