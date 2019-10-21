\name{nncross.ppx}
\alias{nncross.ppx}
\title{Nearest Neighbours Between Two Patterns in Any Dimensions}
\description{
  Given two point patterns \code{X} and \code{Y} in
  many dimensional space,
  finds the nearest neighbour in \code{Y} of each point of \code{X}.
}
\usage{
  \method{nncross}{ppx}(X, Y,
          iX=NULL, iY=NULL,
          what = c("dist", "which"),
          \dots,
          k = 1)
}
\arguments{
  \item{X,Y}{
    Point patterns in any number of spatial dimensions
    (objects of class \code{"ppx"}).
  }
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
  \item{\dots}{Ignored.}
}
\details{
  Given two point patterns \code{X} and \code{Y} in
  \eqn{m}-dimensional space,
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
  B <- boxx(c(0,1), c(0,1), c(0,1), c(0,1))
  ## two different point patterns
  X <- runifpointx(5, B)
  Y <- runifpointx(10, B)
  nncross(X,Y)
  N23 <- nncross(X,Y, k=2:3)

  ## two patterns with some points in common
  Z <- runifpointx(20, B)
  X <- Z[1:15]
  Y <- Z[10:20]
  iX <- 1:15
  iY <- 10:20
  N <- nncross(X,Y, iX, iY, what="which")
  N4 <- nncross(X,Y, iX, iY, k=4)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
