\name{crossdist.pp3} 
\alias{crossdist.pp3}
\title{Pairwise distances between two different three-dimensional point patterns}
\description{
  Computes the distances between pairs of points
  taken from two different three-dimensional point patterns.
}
\usage{
  \method{crossdist}{pp3}(X, Y, \dots, periodic=FALSE, squared=FALSE)
}
\arguments{
  \item{X,Y}{
    Point patterns in three dimensions (objects of class \code{"pp3"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{periodic}{
    Logical. Specifies whether to apply a periodic edge correction.
  }
  \item{squared}{
    Logical. If \code{squared=TRUE}, the squared distances are
    returned instead (this computation is faster).
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th point in \code{X}
  to the \code{j}-th point in \code{Y}.
}
\details{
  Given two point patterns in three-dimensional space, 
  this function computes the Euclidean distance from each point
  in the first pattern to each point in the second pattern,
  and returns a matrix containing these distances.

  This is a method for the generic function \code{\link{crossdist}}
  for three-dimensional point patterns (objects of class \code{"pp3"}).

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
}
\seealso{
  \code{\link{crossdist}},
  \code{\link{pairdist}},
  \code{\link{nndist}},
  \code{\link{G3est}}
}
\examples{
   X <- runifpoint3(20)
   Y <- runifpoint3(30)
   d <- crossdist(X, Y)
   d <- crossdist(X, Y, periodic=TRUE)
}
\author{
  \adrian
  
  
  based on code for two dimensions by 
  Pavel Grabarnik.
}
\keyword{spatial}
\keyword{math}

