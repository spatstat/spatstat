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
linnet(vertices, m, edges)
}
\arguments{
  \item{vertices}{
    Point pattern (object of class \code{"ppp"})
    specifying the vertices of the network.
  }
  \item{m}{
    Adjacency matrix. A matrix of logical values equal to \code{TRUE}
    when the corresponding vertices are joined by a line.
    (Specify either \code{m} or \code{edges}.)
  }
  \item{edges}{
    Edge list. A two-column matrix of integers,
    specifying all pairs of vertices
    that should be joined by an edge. 
    (Specify either \code{m} or \code{edges}.)
  }
}
\details{
  An object of class \code{"linnet"} represents a network of
  straight line segments in two dimensions. The function \code{linnet} creates
  such an object from the minimal information: the spatial location
  of each vertex (endpoint, crossing point or meeting point of lines)
  and information about which vertices are joined by an edge.

  This function can take some time to execute, because the algorithm
  computes various properties of the network that are stored in the
  resulting object.
}
\value{
  Object of class \code{"linnet"} representing the linear network.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{simplenet}} for an example of a linear network.

  \code{\link[spatstat:methods.linnet]{methods.linnet}} for
  methods applicable to \code{linnet} objects.
  
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
