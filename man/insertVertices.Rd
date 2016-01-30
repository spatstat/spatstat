\name{insertVertices}
\alias{insertVertices}
\title{
  Insert New Vertices in a Linear Network
}
\description{
  Adds new vertices to a linear network
  at specified locations along the network.
}
\usage{
insertVertices(L, \dots)
}
\arguments{
  \item{L}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{as.lpp}}
    specifying the positions of the new vertices along the network.
  }
}
\details{
  The new vertices must be locations along the existing network
  \code{L}. 

  The required locations can be specified either as a
  point pattern (class \code{"lpp"} or \code{"ppp"})
  or using coordinate vectors \code{x,y} or \code{seg,tp}
  or \code{x,y,seg,tp} as explained in the help for \code{\link{as.lpp}}.

  This function breaks the existing line segments
  of \code{L} into pieces at the locations specified by
  the coordinates \code{seg,tp} and creates new vertices at these
  locations.

  The result is the modified network, with an attribute \code{"id"} such that
  the point with coordinates \code{seg[i], tp[i]} has become the
  \code{id[i]}th vertex of the new network
}
\value{
  A linear network.
}
\author{
  Adrian Baddeley
}
\seealso{
  \code{\link{as.lpp}}
}
\examples{
   opa <- par(mfrow=c(1,3), mar=rep(0,4))
   simplenet

   plot(simplenet, main="")
   plot(vertices(simplenet), add=TRUE)

   # add two new vertices at specified local coordinates
   L <- insertVertices(simplenet, seg=c(3,7), tp=c(0.2, 0.5))
   L
   plot(L, main="")
   plot(vertices(L), add=TRUE)
   id <- attr(L, "id")
   id
   plot(vertices(L)[id], add=TRUE, pch=16)

   # add new vertices at three randomly-generated points
   X <- runiflpp(3, L)
   LL <- insertVertices(simplenet, X)
   plot(LL, main="")
   plot(vertices(LL), add=TRUE)
   ii <- attr(LL, "id")
   plot(vertices(LL)[ii], add=TRUE, pch=16)
   par(opa)
}
\keyword{spatial}
\keyword{manip}

