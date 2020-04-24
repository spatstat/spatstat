\name{nndist.lpp}
\alias{nndist.lpp}
\title{
  Nearest neighbour distances on a linear network
}
\description{
  Given a pattern of points on a linear network, compute the
  nearest-neighbour distances, measured
  by the shortest path in the network.
}
\usage{
\method{nndist}{lpp}(X, ..., k=1, by=NULL, method="C")
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{k}{
    Integer, or integer vector. The algorithm will compute the distance to the
    \code{k}th nearest neighbour. 
  }
  \item{by}{
    Optional. A factor, which separates \code{X} into groups.
    The algorithm will compute the distance to
    the nearest point in each group. 
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
  nearest neighbour distance for each point (i.e. the distance
  from each point to the nearest other point), measuring
  distance by the shortest path in the network.
  
  If \code{method="C"} the distances are computed using
  code in the C language. If \code{method="interpreted"} then the
  computation is performed using interpreted \R code. The \R code is
  much slower, but is provided for checking purposes.

  The \code{k}th nearest neighbour distance is infinite
  if the \code{k}th nearest neighbour does not exist. This can occur
  if there are fewer than \code{k+1} points in the dataset, or if
  the linear network is not connected.

  If the argument \code{by} is given, it should be a \code{factor},
  of length equal to the number of points in \code{X}.
  This factor effectively partitions \code{X} into subsets,
  each subset associated with one of the levels of \code{X}.
  The algorithm will then compute, for each point of \code{X},
  the distance to the nearest neighbour \emph{in each subset}.
}
\value{
  A numeric vector, of length equal to the number of points in \code{X},
  or a matrix, with one row for each point in \code{X} and one column
  for each entry of \code{k}. Entries are nonnegative numbers or
  infinity (\code{Inf}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{lpp}}
}
\examples{
   X <- runiflpp(12, simplenet)
   nndist(X)
   nndist(X, k=2)

   marks(X) <- factor(rep(letters[1:3], 4))
   nndist(X, by=marks(X))
}
\keyword{spatial}
