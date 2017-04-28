\name{lineardirichlet}
\alias{lineardirichlet}
\title{
  Dirichlet Tessellation on a Linear Network
}
\description{
  Given a point pattern on a linear network, compute the Dirichlet
  (or Voronoi or Thiessen) tessellation induced by the points.
}
\usage{
lineardirichlet(X)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
}
\details{
  The Dirichlet tessellation induced by a point pattern \code{X}
  on a linear network \code{L}
  is a partition of \code{L} into subsets. The subset \code{L[i]}
  associated with the data point \code{X[i]} is the part of \code{L}
  lying closer to \code{X[i]} than to any other data point \code{X[j]},
  where distance is measured by the shortest path.
}
\section{Missing tiles}{
  If the linear network is not connected, and if one of the connected
  components contains no data points, then the Dirichlet tessellation
  is mathematically undefined inside this component.
  The resulting tessellation object includes
  a tile with label \code{NA}, which contains this component of the network.
  A plot of the tessellation will not show this tile. 
}
\value{
  A tessellation on a linear network
  (object of class \code{"lintess"}).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{lintess}}
}
\examples{
  X <- runiflpp(5, simplenet)
  plot(lineardirichlet(X), lwd=3)
  points(X)
}
\keyword{spatial}
\keyword{manip}
