\name{crossdist.ppx}
\alias{crossdist.ppx}
\title{Pairwise Distances Between Two Different Multi-Dimensional Point Patterns}
\description{
  Computes the distances between pairs of points
  taken from two different multi-dimensional point patterns.
}
\usage{
  \method{crossdist}{ppx}(X, Y, \dots)
}
\arguments{
  \item{X,Y}{
    Multi-dimensional point patterns (objects of class \code{"ppx"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{coords.ppx}} to determine
    which coordinates should be used.
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th point in \code{X}
  to the \code{j}-th point in \code{Y}.
}
\details{
  Given two point patterns in multi-dimensional space, 
  this function computes the Euclidean distance from each point
  in the first pattern to each point in the second pattern,
  and returns a matrix containing these distances.

  This is a method for the generic function \code{\link{crossdist}}
  for three-dimensional point patterns (objects of class \code{"ppx"}).

  This function expects two multidimensional
  point patterns \code{X} and \code{Y}, and returns the matrix
  whose \code{[i,j]} entry is the distance from \code{X[i]} to
  \code{Y[j]}.
  
  By default, both spatial and temporal coordinates are extracted.
  To obtain the spatial distance between points in a space-time point
  pattern, set \code{temporal=FALSE}.
}
\seealso{
  \code{\link{crossdist}},
  \code{\link{pairdist}},
  \code{\link{nndist}}
}
\examples{
   df <- data.frame(x=runif(3),y=runif(3),z=runif(3),w=runif(3))
   X <- ppx(data=df)
   df <- data.frame(x=runif(5),y=runif(5),z=runif(5),w=runif(5))
   Y <- ppx(data=df)
   d <- crossdist(X, Y)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\keyword{spatial}
\keyword{math}

