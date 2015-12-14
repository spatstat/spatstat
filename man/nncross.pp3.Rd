\name{nncross.pp3}
\alias{nncross.pp3}
\title{Nearest Neighbours Between Two Patterns in 3D}
\description{
  Given two point patterns \code{X} and \code{Y} in three dimensions,
  finds the nearest neighbour in \code{Y} of each point of \code{X}.
}
\usage{
  \method{nncross}{pp3}(X, Y,
          iX=NULL, iY=NULL,
          what = c("dist", "which"),
          \dots,
          k = 1,
          sortby=c("range", "var", "x", "y", "z"),
          is.sorted.X = FALSE,
          is.sorted.Y = FALSE)
}
\arguments{
  \item{X,Y}{Point patterns in three dimensions
    (objects of class \code{"pp3"}).}
  \item{iX, iY}{Optional identifiers, 
    used to determine whether a point in
    \code{X} is identical to a point in \code{Y}. See Details.
  }
  \item{what}{
    Character string specifying what information should be returned.
    Either the nearest neighbour distance (\code{"dist"}),
    the identifier of the nearest neighbour (\code{"which"}),
    or both.
  }
  \item{k}{
    Integer, or integer vector. The algorithm will compute the distance to the
    \code{k}th nearest neighbour.
  }
  \item{sortby}{
    Determines which coordinate to use to sort the point patterns.
    See Details.
  }
  \item{is.sorted.X, is.sorted.Y}{
    Logical values attesting whether the point patterns \code{X} and
    \code{Y} have been sorted. See Details.
  }
  \item{\dots}{Ignored.}
}
\details{
  Given two point patterns \code{X} and \code{Y} in three dimensions,
  this function finds, for each point of \code{X}, 
  the nearest point of \code{Y}. The distance between these points
  is also computed.
  If the argument \code{k} is specified, then the \code{k}-th nearest
  neighbours will be found.

  The return value is a data frame, with rows corresponding to
  the points of \code{X}.  The first column gives the nearest neighbour
  distances (i.e. the \code{i}th entry is the distance 
  from the \code{i}th point of \code{X} to the nearest element of
  \code{Y}). The second column gives the indices of the nearest
  neighbours (i.e.\ the \code{i}th entry is the index of
  the nearest element in \code{Y}.)
  If \code{what="dist"} then only the vector of distances is returned.
  If \code{what="which"} then only the vector of indices is returned.

  The argument \code{k} may be an integer or an integer vector.
  If it is a single integer, then the \code{k}-th nearest neighbours
  are computed. If it is a vector, then the \code{k[i]}-th nearest
  neighbours are computed for each entry \code{k[i]}. For example, setting
  \code{k=1:3} will compute the nearest, second-nearest and
  third-nearest neighbours. The result is a data frame.

  Note that this function is not symmetric in \code{X} and \code{Y}.
  To find the nearest neighbour in \code{X} of each point in \code{Y},
  use \code{nncross(Y,X)}.

  The arguments \code{iX} and \code{iY} are used when
  the two point patterns \code{X} and \code{Y} have some points in
  common.  In this situation \code{nncross(X, Y)} would return some zero
  distances. To avoid this, attach a unique integer identifier to
  each point, such that two points are identical if their
  identifying numbers are equal. Let \code{iX} be the vector of
  identifier values for the points in \code{X}, and \code{iY}
  the vector of identifiers for points in \code{Y}. Then the code
  will only compare two points if they have different values of the
  identifier. See the Examples.
}
\section{Sorting data and pre-sorted data}{
  Read this section if you care about the speed of computation.
  
  For efficiency, the algorithm sorts both
  the point patterns \code{X} and \code{Y}
  into increasing order of the \eqn{x} coordinate,
  or both into increasing order of the \eqn{y} coordinate,
  or both into increasing order of the \eqn{z} coordinate.
  Sorting is only an intermediate step;
  it does not affect the output, which is always given in the same
  order as the original data.
  
  By default (if \code{sortby="range"}),
  the sorting will occur on the coordinate that has the largest range of
  values (according to the frame of the enclosing window of \code{Y}).
  If \code{sortby = "var"}), sorting will occur on the coordinate that
  has the greater variance (in the pattern \code{Y}).
  Setting \code{sortby="x"} or \code{sortby = "y"}
  or \code{sortby = "z"} will specify that
  sorting should occur on the \eqn{x}, \eqn{y} or \eqn{z} coordinate,
  respectively.

  If the point pattern \code{X} is already
  sorted, then the corresponding argument \code{is.sorted.X}
  should be set to \code{TRUE}, and \code{sortby} should be set
  equal to \code{"x"}, \code{"y"} or \code{"z"} to indicate which coordinate
  is sorted.

  Similarly if \code{Y} is already sorted, then \code{is.sorted.Y}
  should be set to \code{TRUE}, and \code{sortby} should be set
  equal to \code{"x"}, \code{"y"} or \code{"z"} to indicate which coordinate
  is sorted.

  If both \code{X} and \code{Y} are sorted \emph{on the same coordinate
  axis} then both \code{is.sorted.X} and \code{is.sorted.Y}
  should be set to \code{TRUE}, and \code{sortby} should be set
  equal to \code{"x"}, \code{"y"} or \code{"z"} to indicate which coordinate
  is sorted.  
}
\value{
  A data frame, or a vector if the data frame would contain only one column.
  
  By default (if \code{what=c("dist", "which")} and \code{k=1})
  a data frame with two columns:
  \item{dist}{Nearest neighbour distance}
  \item{which}{Nearest neighbour index in \code{Y}}

  If \code{what="dist"} and \code{k=1}, a vector of nearest neighbour distances.

  If \code{what="which"} and \code{k=1}, a vector of nearest neighbour
  indices.

  If \code{k} is specified, the result is a data frame with
  columns containing the \code{k}-th nearest neighbour distances
  and/or nearest neighbour indices.
}
\seealso{
  \code{\link{nndist}} for nearest neighbour
  distances in a single point pattern.
}
\examples{
  # two different point patterns
  X <- pp3(runif(10), runif(10), runif(10), box3(c(0,1)))
  Y <- pp3(runif(20), runif(20), runif(20), box3(c(0,1)))
  N <- nncross(X,Y)$which
  N <- nncross(X,Y, what="which") #faster
  # note that length(N) = 10

  # k-nearest neighbours
  N3 <- nncross(X, Y, k=1:3)

  # two patterns with some points in common
  Z <- pp3(runif(20), runif(20), runif(20), box3(c(0,1)))
  X <- Z[1:15]
  Y <- Z[10:20]
  iX <- 1:15
  iY <- 10:20
  N <- nncross(X,Y, iX, iY, what="which")

}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  ,
  Rolf Turner \email{r.turner@auckland.ac.nz}
  ,
  and Jens Oehlschlaegel
}
\keyword{spatial}
\keyword{math}
