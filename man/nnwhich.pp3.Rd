\name{nnwhich.pp3}
\alias{nnwhich.pp3}
\title{Nearest neighbours in three dimensions}
\description{
  Finds the nearest neighbour of each point in a three-dimensional
  point pattern.
}
\usage{
  \method{nnwhich}{pp3}(X, \dots, k=1)
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
  Numeric vector or matrix giving, for each point,
  the index of its nearest neighbour (or \code{k}th nearest neighbour).

  If \code{k = 1} (the default), the return value is a
  numeric vector \code{v} giving the indices of the nearest neighbours
  (the nearest neighbout of the \code{i}th point is
  the \code{j}th point where \code{j = v[i]}).
  
  If \code{k} is a single integer, then the return value is a
  numeric vector giving the indices of the
  \code{k}th nearest neighbours.

  If \code{k} is a vector, then the return value is a
  matrix \code{m} such that \code{m[i,j]} is the
  index of the \code{k[j]}th nearest neighbour for the
  \code{i}th data point.
}
\details{
  For each point in the given three-dimensional
  point pattern, this function finds
  its nearest neighbour (the nearest other point of the pattern).
  By default it returns a vector giving, for each point,
  the index of the point's
  nearest neighbour. If \code{k} is specified, the algorithm finds
  each point's \code{k}th nearest neighbour.

  The function \code{nnwhich} is generic. This is the method
  for the class \code{"pp3"}.

  If there are no points in the pattern,
  a numeric vector of length zero is returned.
  If there is only one point,
  then the nearest neighbour is undefined, and a value of \code{NA}
  is returned. In general if the number of points is less than or equal
  to \code{k}, then a vector of \code{NA}'s is returned.

  To evaluate the \emph{distance} between a point and its nearest
  neighbour, use \code{\link{nndist}}.

  To find the nearest neighbours from one point pattern
  to another point pattern, use \code{\link{nncross}}.
}
\section{Warnings}{
  A value of \code{NA} is returned if there is only one point
  in the point pattern. 
}
\seealso{
  \code{\link{nnwhich}},
  \code{\link{nndist}},
  \code{\link{nncross}}
}
\examples{
   X <- runifpoint3(30)
   m <- nnwhich(X)
   m2 <- nnwhich(X, k=2)
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  based on two-dimensional code by
  Pavel Grabarnik
}
\keyword{spatial}
\keyword{math}

