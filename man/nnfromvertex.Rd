\name{nnfromvertex}
\alias{nnfromvertex}
\title{
  Nearest Data Point From Each Vertex in a Network
}
\description{
  Given a point pattern on a linear network, 
  for each vertex of the network find the nearest data point.
}
\usage{
nnfromvertex(X, what = c("dist", "which"), k = 1)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{what}{
    Character string specifying whether to return
    the nearest-neighbour distances, nearest-neighbour identifiers,
    or both.
  }
  \item{k}{
    Integer, or integer vector, specifying
    that the \code{k}th nearest neighbour should be returned.
  }
}
\details{
  For each vertex (node) of the linear network,
  this algorithm finds the nearest data point to the vertex,
  and returns either the distance from the vertex to its nearest
  neighbour in \code{X}, or the serial number of the nearest neighbour
  in \code{X}, or both.

  If \code{k} is an integer, then the \code{k}-th nearest neighbour
  is found instead.

  If \code{k} is an integer vector, this is repeated for each
  integer in \code{k}.
}
\value{
  A numeric vector, matrix, or data frame.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{nndist.lpp}}
}
\examples{
  X <- runiflpp(5, simplenet)
  nnfromvertex(X)
  nnfromvertex(X, k=1:3)
}
\keyword{spatial}
\keyword{math}

