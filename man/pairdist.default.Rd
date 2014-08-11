\name{pairdist.default}
\alias{pairdist.default}
\title{Pairwise distances}
\description{
  Computes the matrix of distances between all pairs of points
  in a set of points
}
\usage{
  \method{pairdist}{default}(X, Y=NULL, \dots, period=NULL, method="C", squared=FALSE)
}
\arguments{
  \item{X,Y}{
    Arguments specifying the coordinates of a set of points.
    Typically \code{X} and \code{Y} would be
    numeric vectors of equal length.
    Alternatively \code{Y} may be omitted and \code{X} may be
    a list with two components \code{x} and \code{y},
    or a matrix with two columns.
  }
  \item{\dots}{
    Ignored.
  }
  \item{period}{
    Optional. Dimensions for periodic edge correction.
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
  Given the coordinates of a set of points,
  this function computes the Euclidean distances between all pairs of
  points, and returns the matrix of distances.
  It is a method for the generic function \code{pairdist}.

  The arguments \code{X} and \code{Y} must determine
  the coordinates of a set of points. Typically \code{X} and
  \code{Y} would be numeric vectors of equal length. Alternatively
  \code{Y} may be omitted and \code{X} may be a list with two components
  named \code{x} and \code{y}, or a matrix or data frame with two columns.
  
  Alternatively if \code{period} is given,
  then the distances will be computed in the `periodic'
  sense (also known as `torus' distance).
  The points will be treated as if they are in a rectangle
  of width \code{period[1]} and height \code{period[2]}.
  Opposite edges of the rectangle are regarded as equivalent.
  
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
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link{Kest}}
}
\examples{
   x <- runif(100)
   y <- runif(100)
   d <- pairdist(x, y)
   d <- pairdist(cbind(x,y))
   d <- pairdist(x, y, period=c(1,1))
   d <- pairdist(x, y, squared=TRUE)
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

