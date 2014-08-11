\name{crossdist.ppp} %DontDeclareMethods
\alias{crossdist.ppp}
\title{Pairwise distances between two different point patterns}
\description{
  Computes the distances between pairs of points
  taken from two different point patterns.
}
\usage{
  \method{crossdist}{ppp}(X, Y, \dots, periodic=FALSE, method="C")
}
\arguments{
  \item{X,Y}{
    Point patterns (objects of class \code{"ppp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{periodic}{
    Logical. Specifies whether to apply a periodic edge correction.
  }
  \item{method}{String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th point in \code{X}
  to the \code{j}-th point in \code{Y}.
}
\details{
  Given two point patterns, 
  this function computes the Euclidean distance from each point
  in the first pattern to each point in the second pattern,
  and returns a matrix containing these distances.

  This is a method for the generic function \code{\link{crossdist}}
  for point patterns (objects of class \code{"ppp"}).

  This function expects two
  point patterns \code{X} and \code{Y}, and returns the matrix
  whose \code{[i,j]} entry is the distance from \code{X[i]} to
  \code{Y[j]}.

  Alternatively if \code{periodic=TRUE}, then provided the windows
  containing \code{X} and \code{Y} are identical and are rectangular,
  then the distances will be computed in the `periodic'
  sense (also known as `torus' distance): opposite edges of the
  rectangle are regarded as equivalent.
  This is meaningless if the window is not a rectangle.
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. 
  The C code is faster by a factor of 4.
}
\seealso{
  \code{\link{crossdist}},
  \code{\link{crossdist.default}},
  \code{\link{crossdist.psp}},
  \code{\link{pairdist}},
  \code{\link{nndist}},
  \code{\link{Gest}}
}
\examples{
   data(cells)
   d <- crossdist(cells, runifpoint(6))
   d <- crossdist(cells, runifpoint(6), periodic=TRUE)
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

