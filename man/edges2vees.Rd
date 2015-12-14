\name{edges2vees}
\alias{edges2vees}
\title{
  List Dihedral Triples in a Graph
}
\description{
  Given a list of edges between vertices,
  compile a list of all \sQuote{vees} or dihedral triples
  formed by these edges.
}
\usage{
edges2vees(iedge, jedge, nvert=max(iedge, jedge), \dots,
                check=TRUE)
}
\arguments{
  \item{iedge,jedge}{
    Integer vectors, of equal length, specifying the edges.
  }
  \item{nvert}{
    Number of vertices in the network.
  }
  \item{\dots}{Ignored}
  \item{check}{Logical. Whether to check validity of input data.}
}
\details{
  Given a finite graph with \code{nvert} vertices and with edges
  specified by \code{iedge, jedge}, this low-level function
  finds all \sQuote{vees} or \sQuote{dihedral triples}
  in the graph, that is, all triples
  of vertices \code{(i,j,k)} where \code{i} and \code{j} are joined by
  an edge and \code{i} and \code{k} are joined by an edge. 

  The interpretation of \code{iedge, jedge} is that each successive
  pair of entries specifies an edge in the graph.
  The \eqn{k}th edge joins vertex \code{iedge[k]} to vertex \code{jedge[k]}.
  Entries of \code{iedge} and \code{jedge} must be integers
  from 1 to \code{nvert}.
}
\value{
  A 3-column matrix of integers, in which each row represents a triple
  of vertices, with the first vertex joined to the other two vertices.
}
\seealso{
  \code{\link{edges2triangles}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\examples{
 i <- c(1, 2, 5, 5, 1, 4, 2)
 j <- c(2, 3, 3, 1, 3, 2, 5)
 edges2vees(i, j)
}
\keyword{spatial}
\keyword{manip}
