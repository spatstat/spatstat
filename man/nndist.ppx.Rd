\name{nndist.ppx}
\alias{nndist.ppx}
\title{Nearest Neighbour Distances in Any Dimensions}
\description{
  Computes the distance from each point to its nearest neighbour
  in a multi-dimensional point pattern.
  Alternatively computes the distance to the
  second nearest neighbour, or third nearest, etc.
}
\usage{
  \method{nndist}{ppx}(X, \dots, k=1)
}
\arguments{
  \item{X}{
    Multi-dimensional point pattern
    (object of class \code{"ppx"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{coords.ppx}} to determine
    which coordinates should be used.
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
  in a multi-dimensional
  point pattern to its nearest neighbour (the nearest other
  point of the pattern). If \code{k} is specified, it computes the
  distance to the \code{k}th nearest neighbour.

  The function \code{nndist} is generic; this function
  \code{nndist.ppx} is the method for the class \code{"ppx"}.

  The argument \code{k} may be a single integer, or an integer vector.
  If it is a vector, then the \eqn{k}th nearest neighbour distances are
  computed for each value of \eqn{k} specified in the vector.

  If there is only one point (if \code{x} has length 1),
  then a nearest neighbour distance of \code{Inf} is returned.
  If there are no points (if \code{x} has length zero)
  a numeric vector of length zero is returned.

  To identify \emph{which} point is the nearest neighbour of a given point,
  use \code{\link{nnwhich}}.

  To find the nearest neighbour distances from one point pattern
  to another point pattern, use \code{\link{nncross}}.
  
  By default, both spatial and temporal coordinates are extracted.
  To obtain the spatial distance between points in a space-time point
  pattern, set \code{temporal=FALSE}.
}
\section{Warnings}{
  An infinite or \code{NA} value is returned if the
  distance is not defined (e.g. if there is only one point
  in the point pattern).
}
\seealso{
  \code{\link{nndist}},
  \code{\link{pairdist}},
  \code{\link{nnwhich}}
}
\examples{
   df <- data.frame(x=runif(5),y=runif(5),z=runif(5),w=runif(5))
   X <- ppx(data=df)

   # nearest neighbours
   d <- nndist(X)

   # second nearest neighbours
   d2 <- nndist(X, k=2)

   # first, second and third nearest
   d1to3 <- nndist(X, k=1:3)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\keyword{spatial}
\keyword{math}

