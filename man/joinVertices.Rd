\name{joinVertices}
\alias{joinVertices}
\title{
  Join Vertices in a Network
}
\description{
  Join the specified vertices in a linear network, creating a new network.
}
\usage{
joinVertices(L, from, to)
}
\arguments{
  \item{L}{
    A linear network (object of class \code{"linnet"})
    or point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{from,to}{
    Integers, or integer vectors of equal length,
    specifying the vertices which should be joined.
    Alternatively \code{from} can be a 2-column matrix of integers
    and \code{to} is missing or \code{NULL}.
  }
}
\details{
  Vertices of the network are numbered by their order of appearance
  in the point pattern \code{vertices(L)}.

  If \code{from} and \code{to} are single integers, then the
  pair of vertices numbered \code{from} and \code{to} will be
  joined to make a new segment of the network.
  If \code{from} and \code{to} are vectors of integers, then
  vertex \code{from[i]} will be joined to vertex \code{to[i]} for
  each \code{i = 1,2,..}.

  If \code{L} is a network (class \code{"linnet"}), the result is
  another network, created by adding new segments.
  If \code{L} is a point pattern on a network (class \code{"lpp"}), the
  result is another point pattern object, created by adding new segments to
  the underlying network, and retaining the points.

  In the resulting object, the new line segments are appended to the
  existing list of line segments.
}
\value{
  A linear network (object of class \code{"linnet"})
  or point pattern on a linear network (object of class \code{"lpp"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{linnet}}, \code{\link{methods.linnet}},
  \code{\link{thinNetwork}}
}
\examples{
  snet <- joinVertices(simplenet, 4, 8)
  plot(solist(simplenet, snet), main="")
  X <- runiflpp(3, simplenet)
  Y <- joinVertices(X, 4, 8)
}
\keyword{spatial}
\keyword{manip}
