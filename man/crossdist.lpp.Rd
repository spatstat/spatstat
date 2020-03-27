\name{crossdist.lpp} 
\alias{crossdist.lpp}
\title{Pairwise distances between two point patterns on a linear network}
\description{
  Computes the distances between pairs of points
  taken from two different point patterns on
  the same linear network.
}
\usage{
  \method{crossdist}{lpp}(X, Y, \dots, method="C", check=TRUE)
}
\arguments{
  \item{X,Y}{
    Point patterns on a linear network (objects of class \code{"lpp"}).
    They must lie on the \emph{same} network.
  }
  \item{\dots}{
    Ignored.
  }
  \item{method}{
    String specifying which method of calculation to use
    when the network data use the non-sparse representation.
    Values are \code{"C"} and \code{"interpreted"}.
  }
  \item{check}{
    Logical value specifying whether to check that \code{X} and \code{Y}
    are defined on the same network. Default is \code{check=TRUE}.
    Setting \code{check=FALSE} will save time,
    but should only be used if it is certain that the
    two networks are identical.
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th point in \code{X}
  to the \code{j}-th point in \code{Y}.
  Matrix entries are nonnegative numbers or infinity (\code{Inf}).
}
\details{
  Given two point patterns on a linear network,
  this function computes the distance from each point
  in the first pattern to each point in the second pattern,
  measuring distance by the shortest path along the network.

  This is a method for the generic function \code{\link{crossdist}}
  for the class of point patterns on a linear network
  (objects of class \code{"lpp"}).

  This function expects two point pattern objects \code{X} and \code{Y}
  on the \emph{same} linear network,
  and returns the matrix whose \code{[i,j]} entry is the
  shortest-path distance from \code{X[i]} to \code{Y[j]}.

  If two points cannot be joined by a path,
  the distance between them is infinite (\code{Inf}).

  The argument \code{method} is not normally used. It is
  retained only for developers to check the validity of the software.
}
\section{Algorithms and accuracy}{
  Distances are accurate within the numerical tolerance of the
  network, \code{summary(X)$toler}.

  For network data stored in the non-sparse
  representation described in \code{\link{linnet}},
  then pairwise distances are computed using the matrix of path distances
  between vertices of the network, using \R code if
  \code{method = "interpreted"}, or using \proglang{C} code if 
  \code{method="C"} (the default).
  
  For networks stored in the sparse representation,
  the argument \code{method} has no effect, and the distances are
  computed using an efficient \proglang{C} algorithm.
}
\seealso{
  \code{\link{crossdist}},
  \code{\link{crossdist.ppp}},
  \code{\link{pairdist}},
  \code{\link{nndist}}
}
\examples{
   v <- split(chicago)
   X <- v$cartheft
   Y <- v$burglary
   d <- crossdist(X, Y)
   d[1:3,1:4]
}
\author{
  \adrian.
}
\keyword{spatial}
\keyword{math}

