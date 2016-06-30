\name{nncross.lpp}
\alias{nncross.lpp}
\title{Nearest Neighbours on a Linear Network}
\description{
  Given two point patterns \code{X} and \code{Y} on a linear network,
  finds the nearest neighbour in \code{Y} of each point of \code{X}
  using the shortest path in the network.
}
\usage{
  \method{nncross}{lpp}(X, Y,
          iX=NULL, iY=NULL,
          what = c("dist", "which"),
          \dots,
          k = 1,
          method="C")
}
\arguments{
  \item{X,Y}{
    Point patterns on a linear network (objects of class \code{"lpp"}).
    They must lie on the \emph{same} linear network.
  }
  \item{iX, iY}{
    Optional identifiers, used to determine whether a point in
    \code{X} is identical to a point in \code{Y}. See Details.
  }
  \item{what}{
    Character string specifying what information should be returned.
    Either the nearest neighbour distance (\code{"dist"}),
    the identifier of the nearest neighbour (\code{"which"}),
    or both.
  }
  \item{\dots}{Ignored.}
  \item{k}{
    Integer, or integer vector. The algorithm will compute the distance to the
    \code{k}th nearest neighbour, for each value of \code{k}.
  }
  \item{method}{
    Internal use only.
  }
}
\details{
  Given two point patterns \code{X} and \code{Y} on the same linear
  network, this function finds, for each point of \code{X}, 
  the nearest point of \code{Y}, measuring distance by the shortest path
  in the network. The distance between these points
  is also computed.

  The return value is a data frame, with rows corresponding to
  the points of \code{X}.  The first column gives the nearest neighbour
  distances (i.e. the \code{i}th entry is the distance 
  from the \code{i}th point of \code{X} to the nearest element of
  \code{Y}). The second column gives the indices of the nearest
  neighbours (i.e.\ the \code{i}th entry is the index of
  the nearest element in \code{Y}.)
  If \code{what="dist"} then only the vector of distances is returned.
  If \code{what="which"} then only the vector of indices is returned.

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
  By default (if \code{what=c("dist", "which")} and \code{k=1})
  a data frame with two columns:
  \item{dist}{Nearest neighbour distance}
  \item{which}{Nearest neighbour index in \code{Y}}

  If \code{what="dist"}, a vector of nearest neighbour distances.

  If \code{what="which"}, a vector of nearest neighbour indices.

  If \code{k} is a vector of integers, the result is a matrix
  with one row for each point in \code{X},
  giving the distances and/or indices of the \code{k}th nearest
  neighbours in \code{Y}.
}
\seealso{
  \code{\link{nndist.lpp}} for nearest neighbour
  distances in a single point pattern.

  \code{\link{nnwhich.lpp}} to identify which points are nearest
  neighbours in a single point pattern.
}
\examples{
  # two different point patterns
  X <- runiflpp(3, simplenet)
  Y <- runiflpp(5, simplenet)
  nn <- nncross(X,Y)
  nn
  plot(simplenet, main="nncross")
  plot(X, add=TRUE, cols="red")
  plot(Y, add=TRUE, cols="blue", pch=16)
  XX <- as.ppp(X)
  YY <- as.ppp(Y)
  i <- nn$which
  arrows(XX$x, XX$y, YY[i]$x, YY[i]$y, length=0.15)

  # nearest and second-nearest neighbours
  nncross(X, Y, k=1:2)

  # two patterns with some points in common
  X <- Y[1:2]
  iX <- 1:2
  iY <- 1:5
  nncross(X,Y, iX, iY)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{math}
