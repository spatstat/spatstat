\name{edges2triangles}
\alias{edges2triangles}
\title{
  List Triangles in a Graph
}
\description{
  Given a list of edges between vertices,
  compile a list of all triangles formed by these edges.
}
\usage{
edges2triangles(iedge, jedge, nvert=max(iedge, jedge), \dots,
                check=TRUE, friendly=rep(TRUE, nvert))
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
  \item{friendly}{
    Optional. For advanced use. See Details.
  }
}
\details{
  This low level function finds all the triangles (cliques of size 3)
  in a finite graph with \code{nvert} vertices and with edges
  specified by \code{iedge, jedge}. 

  The interpretation of \code{iedge, jedge} is that each successive
  pair of entries specifies an edge in the graph.
  The \eqn{k}th edge joins vertex \code{iedge[k]} to vertex \code{jedge[k]}.
  Entries of \code{iedge} and \code{jedge} must be integers
  from 1 to \code{nvert}.

  To improve efficiency in some applications, the optional argument
  \code{friendly} can be used. It should be a logical vector of
  length \code{nvert} specifying a labelling of the vertices,
  such that two vertices \code{j,k} which are \emph{not} friendly
  (\code{friendly[j] = friendly[k] = FALSE})
  are \emph{never} connected by an edge. 
}
\value{
  A 3-column matrix of integers, in which each row represents a triangle.
}
\seealso{
  \code{\link{edges2vees}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\examples{
 i <- c(1, 2, 5, 5, 1, 4, 2)
 j <- c(2, 3, 3, 1, 3, 2, 5)
 edges2triangles(i, j)
}
\keyword{spatial}
\keyword{manip}
