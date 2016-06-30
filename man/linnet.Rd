\name{linnet}
\alias{linnet}
\title{
  Create a Linear Network
}
\description{
  Creates an object of class \code{"linnet"} representing
  a network of line segments.
}
\usage{
linnet(vertices, m, edges, sparse=FALSE, warn=TRUE)
}
\arguments{
  \item{vertices}{
    Point pattern (object of class \code{"ppp"})
    specifying the vertices of the network.
  }
  \item{m}{
    Adjacency matrix. A matrix or sparse matrix
    of logical values equal to \code{TRUE}
    when the corresponding vertices are joined by a line.
    (Specify either \code{m} or \code{edges}.)
  }
  \item{edges}{
    Edge list. A two-column matrix of integers,
    specifying all pairs of vertices
    that should be joined by an edge. 
    (Specify either \code{m} or \code{edges}.)
  }
  \item{sparse}{
    Optional. Logical value indicating whether to use a
    sparse matrix representation of the network. See Details.
  }
  \item{warn}{
    Logical value indicating whether to issue a warning if the resulting
    network is not connected.
  }
}
\details{
  An object of class \code{"linnet"} represents a network of
  straight line segments in two dimensions. The function \code{linnet} creates
  such an object from the minimal information: the spatial location
  of each vertex (endpoint, crossing point or meeting point of lines)
  and information about which vertices are joined by an edge.

  If \code{sparse=FALSE} (the default), the algorithm will compute
  and store various properties of the network, including
  the adjacency matrix \code{m} and a matrix giving the
  shortest-path distances between each pair of vertices in the network.
  This is more efficient for small datasets. However it can require
  large amounts of memory and can take a long time to execute.

  If \code{sparse=TRUE}, then the shortest-path distances will not be computed,
  and the network adjacency matrix \code{m} will be stored as a
  sparse matrix. This saves a lot of time and memory when creating the
  linear network.

  If the argument \code{edges} is given, then it will also determine
  the \emph{ordering} of the line segments when they are stored or extracted.
  For example, \code{edges[i,]} corresponds to \code{as.psp(L)[i]}.
}
\value{
  Object of class \code{"linnet"} representing the linear network.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  \adrian 
}
\seealso{
  \code{\link{simplenet}} for an example of a linear network.

  \code{\link[spatstat:methods.linnet]{methods.linnet}} for
  methods applicable to \code{linnet} objects.

  Special tools: \code{\link{thinNetwork}}, \code{\link{insertVertices}},
  \code{\link{connected.linnet}}, \code{\link{lixellate}}.
  
  \code{\link{delaunayNetwork}} for the Delaunay triangulation
  as a network.
  
  \code{\link{ppp}},
  \code{\link{psp}}.
}
\examples{
  # letter 'A' specified by adjacency matrix
  v <- ppp(x=(-2):2, y=3*c(0,1,2,1,0), c(-3,3), c(-1,7))
  m <- matrix(FALSE, 5,5)
  for(i in 1:4) m[i,i+1] <- TRUE
  m[2,4] <- TRUE
  m <- m | t(m)
  letterA <- linnet(v, m)
  plot(letterA)

  # letter 'A' specified by edge list
  edg <- cbind(1:4, 2:5)
  edg <- rbind(edg, c(2,4))
  letterA <- linnet(v, edges=edg)
}
\keyword{spatial}
