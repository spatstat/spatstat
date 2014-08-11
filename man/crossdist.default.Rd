\name{crossdist.default}
\alias{crossdist.default}
\title{Pairwise distances between two different sets of points}
\description{
  Computes the distances between each pair of points
  taken from two different sets of points.
}
\usage{
  \method{crossdist}{default}(X, Y, x2, y2, \dots, period=NULL, method="C")
}
\arguments{
  \item{X,Y}{
    Numeric vectors of equal length specifying the coordinates of
    the first set of points.
  }
  \item{x2,y2}{
    Numeric vectors of equal length specifying the coordinates of
    the second set of points.
  }
  \item{\dots}{
    Ignored.
  }
  \item{period}{
    Optional. Dimensions for periodic edge correction.
  }
  \item{method}{String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th point in the first set of points
  to the \code{j}-th point in the second set of points.
}
\details{
  Given two sets of points,
  this function computes the Euclidean distance from each point
  in the first set to each point in the second set,
  and returns a matrix containing these distances.

  This is a method for the generic function \code{\link{crossdist}}.

  This function expects \code{X} and \code{Y} to be numeric vectors
  of equal length specifying the coordinates of the first set of points.
  The arguments \code{x2},\code{y2} specify the coordinates of the
  second set of points.
  
  Alternatively if \code{period} is given,
  then the distances will be computed in the `periodic'
  sense (also known as `torus' distance).
  The points will be treated as if they are in a rectangle
  of width \code{period[1]} and height \code{period[2]}.
  Opposite edges of the rectangle are regarded as equivalent.
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. 
  The C code is faster by a factor of 4.
}
\seealso{
  \code{\link{crossdist}},
  \code{\link{crossdist.ppp}},
  \code{\link{crossdist.psp}},
  \code{\link{pairdist}},
  \code{\link{nndist}},
  \code{\link{Gest}}
}
\examples{
   d <- crossdist(runif(7), runif(7), runif(12), runif(12))
   d <- crossdist(runif(7), runif(7), runif(12), runif(12), period=c(1,1))
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

