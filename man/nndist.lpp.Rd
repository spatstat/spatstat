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
\method{nndist}{lpp}(X, ..., k=1, method="C")
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{method}{
    Optional string determining the method of calculation.
    Either \code{"interpreted"} or \code{"C"}.
  }
  \item{k}{
    Integer, or integer vector. The algorithm will compute the distance to the
    \code{k}th nearest neighbour. 
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
}
\value{
  A numeric vector, of length equal to the number of points in \code{X},
  or a matrix, with one row for each point in \code{X} and one column
  for each entry of \code{k}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
}
\seealso{
  \code{\link{lpp}}
}
\examples{
   X <- runiflpp(12, simplenet)
   nndist(X)
   nndist(X, k=2)
}
\keyword{spatial}
