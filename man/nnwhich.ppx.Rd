\name{nnwhich.ppx}
\alias{nnwhich.ppx}
\title{Nearest Neighbours in Any Dimensions}
\description{
  Finds the nearest neighbour of each point in a multi-dimensional
  point pattern.
}
\usage{
  \method{nnwhich}{ppx}(X, \dots, k=1)
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
  For each point in the given multi-dimensional
  point pattern, this function finds
  its nearest neighbour (the nearest other point of the pattern).
  By default it returns a vector giving, for each point,
  the index of the point's
  nearest neighbour. If \code{k} is specified, the algorithm finds
  each point's \code{k}th nearest neighbour.

  The function \code{nnwhich} is generic. This is the method
  for the class \code{"ppx"}.

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
  
  By default, both spatial and temporal coordinates are extracted.
  To obtain the spatial distance between points in a space-time point
  pattern, set \code{temporal=FALSE}.
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
   df <- data.frame(x=runif(5),y=runif(5),z=runif(5),w=runif(5))
   X <- ppx(data=df)
   m <- nnwhich(X)
   m2 <- nnwhich(X, k=2)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\keyword{spatial}
\keyword{math}

