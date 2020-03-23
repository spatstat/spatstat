\name{addVertices}
\alias{addVertices}
\title{
  Add New Vertices to a Linear Network
}
\description{
  Adds new vertices to a linear network
  at specified locations outside the network.
}
\usage{
addVertices(L, X, join=NULL, joinmarks=NULL)
}
\arguments{
  \item{L}{
    Existing linear network (object of class \code{"linnet"})
    or point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{X}{
    Point pattern (object of class \code{"ppp"})
    specifying the new vertices.
  }
  \item{join}{
    Optional information specifying how to join the new vertices
    \code{X} to the existing network. See Details.
    If \code{join=NULL} (the default),
    the new vertices are simply added to the list of network vertices
    without being joined to the rest of the network.
  }
  \item{joinmarks}{
    Optional vector or data frame of marks associated with the
    new edges specified by \code{join}.
  }
}
\details{
  This function adds new vertices to an existing
  linear network \code{L}, at specified locations \code{X} outside the network.
  
  The argument \code{L} can be either a linear network (class
  \code{"linnet"}) or some other object that includes a linear network.
 
  The new vertex locations are points outside the network,
  specified as a point pattern \code{X} (object of class \code{"ppp"}).

  The argument \code{join} specifies how to join the new vertices
  to the existing network. 
  \itemize{
    \item
    If \code{join=NULL} (the default),
    the new vertices are simply added to the list of network vertices
    without being joined to the rest of the network.
    \item
    If \code{join} is a vector of integers, then these are taken to be
    indices of existing vertices of \code{L} in the order given
    in \code{V = vertices(L)}. Then each new vertex \code{X[i]} will be
    joined to an existing vertex \code{V[j]} where \code{j = join[i]}.
    Each new vertex is joined to exactly one existing vertex.
    \item
    If \code{join="vertices"} then each new vertex \code{X[i]} is joined
    to the nearest existing vertex \code{V[j]}.
    Each new vertex is joined to exactly one existing vertex.
    \item
    If \code{join="nearest"} then each new vertex is projected to the
    nearest location along on the network; these locations are inserted
    as new vertices of \code{L}; and then each vertex \code{X[i]}
    is joined to the corresponding projected point.
    Each new vertex is joined to exactly one newly-inserted vertex.
    \item
    If \code{join} is a point pattern on a network (class \code{"lpp"}),
    it must be defined on the same network as \code{L} and it must
    consist of the same number of points as \code{X}. The points of
    \code{join} will be inserted as new vertices of \code{L},
    and then each vertex \code{X[i]} is joined to the corresponding
    point \code{join[i]}.
    Each new vertex is joined to exactly one newly-inserted vertex.
  }
  The result is the modified object, with an attribute \code{"id"} such that
  the \code{i}th added vertex has become the
  \code{id[i]}th vertex of the new network.
}
\value{
  An object of the same class as \code{L} representing the result of
  adding the new vertices.
  The result also has an attribute \code{"id"} as described in Details.
}
\author{
  Adrian Baddeley
}
\seealso{
  \code{\link{insertVertices}} to insert vertices along an existing network.
  
  \code{\link{as.lpp}},
  \code{\link{linnet}},
  \code{\link{methods.linnet}},
  \code{\link{joinVertices}},
  \code{\link{thinNetwork}}.
}
\examples{
   opa <- par(mfrow=c(1,3))
   L <- simplenet
   X <- runifpoint(20, Window(simplenet))
   plot(L)
   plot(X, add=TRUE, cols="green", pch=16, cex=2)
   plot(addVertices(L, X, "nearest"), col="red")
   plot(L, add=TRUE, col="grey", lwd=3)
   plot(X, add=TRUE, cols="green", pch=16, cex=2)
   plot(addVertices(L, X, "vertices"), col="red")
   plot(L, add=TRUE, col="grey", lwd=3)
   plot(X, add=TRUE, cols="green", pch=16, cex=2)
   par(opa)
}
\keyword{spatial}
\keyword{manip}

