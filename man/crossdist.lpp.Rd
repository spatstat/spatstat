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
  retained only for checking the validity of the software.
  It applies only when the network data are stored in the non-sparse
  representation described in \code{\link{linnet}}.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code. If \code{method="C"}
  (the default) then the distances are computed using C code.
  For networks stored in the sparse representation,
  the argument \code{method} has no effect, and the distances are
  computed using an efficient C algorithm.
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
}
\author{
  \adrian.
}
\keyword{spatial}
\keyword{math}

