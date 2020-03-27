\name{pairdist.lpp}
\alias{pairdist.lpp}
\title{
  Pairwise shortest-path distances between points on
  a linear network
}
\description{
  Given a pattern of points on a linear network, compute the
  matrix of distances between all pairs of points, measuring
  distance by the shortest path in the network.
}
\usage{
\method{pairdist}{lpp}(X, ..., method="C")
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{method}{
    Optional string determining the method of calculation.
    Either \code{"interpreted"} or \code{"C"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  Given a pattern of points on a linear network, this function computes the
  matrix of distances between all pairs of points, measuring
  distance by the shortest path in the network. 
  
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
\value{
  A symmetric matrix, whose values are nonnegative numbers or infinity
  (\code{Inf}).
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  \adrian.
}
\seealso{
  \code{\link{lpp}}
}
\examples{
   X <- runiflpp(12, simplenet)
   d <- pairdist(X)
   d[1:3, 1:3]
}
\keyword{spatial}
