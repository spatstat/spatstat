\name{connected.lpp}
\alias{connected.lpp}
\title{
  Connected Components of a Point Pattern on a Linear Network
}
\description{
  Finds the topologically-connected components of a point pattern on a
  linear network, when all pairs of points closer than a threshold distance
  are joined.
}
\usage{
\method{connected}{lpp}(X, R=Inf, \dots, dismantle=TRUE)
}
\arguments{
  \item{X}{
    A linear network (object of class \code{"lpp"}).
  }
   \item{R}{
    Threshold distance. Pairs of points will be joined together
    if they are closer than \code{R} units apart, measured
    by the shortest path in the network.
    The default \code{R=Inf} implies that points
    will be joined together if they are mutually connected by any
    path in the network.
  }
  \item{dismantle}{
    Logical. If \code{TRUE} (the default), the network itself will be
    divided into its path-connected components using
    \code{\link{connected.linnet}}.
  }
 \item{\dots}{
    Ignored.
  }
}
\details{
  The function \code{connected} is generic. This is the method for
  point patterns on a linear network (objects of class \code{"lpp"}).
  It divides the point pattern \code{X} into one or more groups of points.

  If \code{R=Inf} (the default), then \code{X} is divided into groups
  such that any pair of points in the same group
  can be joined by a path in the network.

  If \code{R} is a finite number, then two points of \code{X} are
  declared to be \emph{R-close} if they lie closer than
  \code{R} units apart, measured by the length of the shortest path in the
  network. Two points are \emph{R-connected} if they 
  can be reached by a series of steps between R-close pairs of
  points of \code{X}. Then \code{X} is divided into groups such that
  any pair of points in the same group is R-connected.

  If \code{dismantle=TRUE} (the default) the algorithm first checks
  whether the network is connected (i.e. whether any pair of vertices
  can be joined by a path in the network), and if not, the network is
  decomposed into its connected components.
}
\value{
  A point pattern (of class \code{"lpp"}) with marks indicating the
  grouping, or a list of such point patterns.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{thinNetwork}}
}
\examples{
   # remove some edges from a network to make it disconnected
   plot(simplenet, col="grey", main="", lty=2)
   A <- thinNetwork(simplenet, retainedges=-c(3,5))
   plot(A, add=TRUE, lwd=2)
   X <- runiflpp(10, A)
   # find the connected components
   cX <- connected(X)
   plot(cX[[1]], add=TRUE, col="blue", lwd=2)
}
\keyword{spatial}
\keyword{manip}
