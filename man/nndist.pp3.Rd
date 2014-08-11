\name{nndist.pp3}
\alias{nndist.pp3}
\title{Nearest neighbour distances in three dimensions}
\description{
  Computes the distance from each point to its nearest neighbour
  in a three-dimensional point pattern.
  Alternatively computes the distance to the
  second nearest neighbour, or third nearest, etc.
}
\usage{
  \method{nndist}{pp3}(X, \dots, k=1)
}
\arguments{
  \item{X}{
    Three-dimensional point pattern
    (object of class \code{"pp3"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{k}{
    Integer, or integer vector. The algorithm will compute the distance to the
    \code{k}th nearest neighbour. 
  }
}
\value{
  Numeric vector or matrix containing the 
  nearest neighbour distances for each point.

  If \code{k = 1} (the default), the return value is a
  numeric vector \code{v} such that \code{v[i]} is the
  nearest neighbour distance for the \code{i}th data point.
  
  If \code{k} is a single integer, then the return value is a
  numeric vector \code{v} such that \code{v[i]} is the
  \code{k}th nearest neighbour distance for the
  \code{i}th data point.

  If \code{k} is a vector, then the return value is a
  matrix \code{m} such that \code{m[i,j]} is the
  \code{k[j]}th nearest neighbour distance for the
  \code{i}th data point.
}
\details{
  This function computes the Euclidean distance from each point
  in a three-dimensional
  point pattern to its nearest neighbour (the nearest other
  point of the pattern). If \code{k} is specified, it computes the
  distance to the \code{k}th nearest neighbour.

  The function \code{nndist} is generic; this function
  \code{nndist.pp3} is the method for the class \code{"pp3"}.

  The argument \code{k} may be a single integer, or an integer vector.
  If it is a vector, then the \eqn{k}th nearest neighbour distances are
  computed for each value of \eqn{k} specified in the vector.

  If there is only one point (if \code{x} has length 1),
  then a nearest neighbour distance of \code{Inf} is returned.
  If there are no points (if \code{x} has length zero)
  a numeric vector of length zero is returned.

  To identify \emph{which} point is the nearest neighbour of a given point,
  use \code{\link{nnwhich}}.

  To use the nearest neighbour distances for statistical inference,
  it is often advisable to use the edge-corrected empirical distribution,
  computed by \code{\link{G3est}}.

  To find the nearest neighbour distances from one point pattern
  to another point pattern, use \code{\link{nncross}}.
}
\section{Warnings}{
  An infinite or \code{NA} value is returned if the
  distance is not defined (e.g. if there is only one point
  in the point pattern).
}
\seealso{
  \code{\link{nndist}},
  \code{\link{pairdist}},
  \code{\link{G3est}},
  \code{\link{nnwhich}}
}
\examples{
   X <- runifpoint3(40)

   # nearest neighbours
   d <- nndist(X)

   # second nearest neighbours
   d2 <- nndist(X, k=2)

   # first, second and third nearest
   d1to3 <- nndist(X, k=1:3)
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  based on code for two dimensions by
  Pavel Grabarnik
}
\keyword{spatial}
\keyword{math}

