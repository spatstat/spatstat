\name{pairdist.ppx}
\alias{pairdist.ppx}
\title{Pairwise Distances in Any Dimensions}
\description{
  Computes the matrix of distances between all pairs of points
  in a multi-dimensional point pattern.
}
\usage{
  \method{pairdist}{ppx}(X, \dots)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppx"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{coords.ppx}} to determine
    which coordinates should be used.
  }
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance
  between the points numbered \code{i} and \code{j}.
}
\details{
  This is a method for the generic function \code{pairdist}.

  Given a multi-dimensional point pattern \code{X}
  (an object of class \code{"ppx"}),
  this function computes the Euclidean distances between all pairs of
  points in \code{X}, and returns the matrix of distances.

  By default, both spatial and temporal coordinates are extracted.
  To obtain the spatial distance between points in a space-time point
  pattern, set \code{temporal=FALSE}.
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{crossdist}},
  \code{\link{nndist}}
}
\examples{
   df <- data.frame(x=runif(4),y=runif(4),z=runif(4),w=runif(4))
   X <- ppx(data=df)
   pairdist(X)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\keyword{spatial}
\keyword{math}

