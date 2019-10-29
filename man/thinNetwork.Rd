\name{thinNetwork}
\alias{thinNetwork}
\title{
  Remove Vertices or Segments from a Linear Network
}
\description{
  Delete some vertices and/or segments from a linear network
  or related object.
}
\usage{
thinNetwork(X, retainvertices, retainedges)
}
\arguments{
  \item{X}{
    A linear network (object of class \code{"linnet"}),
    or a point pattern on a linear network (object of class
    \code{"lpp"}).
  }
  \item{retainvertices}{
    Optional. Subset index specifying which vertices should be retained
    (not deleted). 
  }
  \item{retainedges}{
    Optional. Subset index specifying which edges (segments) should be retained
    (not deleted). 
  }
}
\details{
  This function deletes some of the vertices and edges (segments) in the
  linear network.

  The arguments \code{retainvertices} and \code{retainedges} can be
  any kind of subset index: a vector of positive integers specifying which
  vertices/edges should be retained; a vector of negative integers
  specifying which vertices/edges should be deleted; or a logical vector
  specifying whether each vertex/edge should be retained (\code{TRUE})
  or deleted (\code{FALSE}).

  Vertices are indexed in the same sequence as in
  \code{vertices(as.linnet(X))}.
  Segments are indexed in the same sequence as in
  \code{as.psp(as.linnet(X))}.

  The argument \code{retainedges} has higher precedence than
  \code{retainvertices} in the sense that:
  \itemize{
    \item If \code{retainedges} is given, then
    any vertex which is an endpoint of a retained edge will also be
    retained.
    \item
    If \code{retainvertices} is given and \code{retainedges} is \bold{missing},
    then any segment joining two retained vertices will also be retained.
    \item
    Thus, when both \code{retainvertices} and \code{retainedges} are
    given, it is possible that more vertices will be retained than those
    specified by \code{retainvertices}.
  }

  After the network has been altered, other consequential changes will
  occur, including renumbering of the segments and vertices.
  If \code{X} is a point pattern on a linear network, then data points
  will be deleted if they lie on a deleted edge.
}
\value{
  An object of the same kind as \code{X}.
}
\author{
  \adrian
  and Suman Rakshit.
}
\seealso{
  \code{\link{linnet}} to make a network;

  \code{\link{connected.linnet}} to extract connected components.

  \code{\link{repairNetwork}}.
}
\examples{
   L <- simplenet
   plot(L, main="thinNetwork(L, retainedges=c(-3, -5))")
   text(midpoints.psp(as.psp(L)), labels=1:nsegments(L), pos=3)
   Lsub <- thinNetwork(L, retainedges=c(-3, -5))
   plot(Lsub, add=TRUE, col="blue", lwd=2)
}
\keyword{spatial}
\keyword{manip}
