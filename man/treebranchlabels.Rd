\name{treebranchlabels}
\alias{treebranchlabels}
\title{
  Label Vertices of a Tree by Branch Membership
}
\description{
  Given a linear network which is a tree (acyclic graph),
  this function assigns a label to each vertex, indicating
  its position in the tree.
}
\usage{
  treebranchlabels(L, root = 1)
}
\arguments{
  \item{L}{
    Linear network (object of class \code{"linnet"}).
    The network must have no loops.
  }
  \item{root}{
    Root of the tree. An integer index identifying
    which point in \code{vertices(L)} is the root of the tree.
  }
}
\details{
  The network \code{L} should be a tree, that is, it must have no loops.

  This function computes a character string label for each vertex
  of the network \code{L}. The vertex identified by \code{root}
  (that is, \code{vertices(L)[root]}) is taken as the root of the tree
  and is given the empty label \code{""}.
  \itemize{
    \item If there are several line
    segments which meet at the root vertex, each of these segments is the
    start of a new branch of the tree; the other endpoints of these
    segments are assigned the labels 
    \code{"a"}, \code{"b"}, \code{"c"} and so on.
    \item If only one segment issues from the root vertex,
    the other endpoint of this segment is assigned the empty label
    \code{""}.
  }
  A similar rule is then applied to each of the newly-labelled vertices.
  If the vertex labelled \code{"a"} is joined to two other unlabelled
  vertices, these will be labelled \code{"aa"} and \code{"ab"}.
  The rule is applied recursively until all vertices have been labelled.

  If \code{L} is not a tree, the algorithm will terminate, but the
  results will be nonsense.
}
\value{
  A vector of character strings, with one entry for each point in
  \code{vertices(L)}.
}
\author{
\spatstatAuthors
}
\seealso{
  \code{\link{deletebranch}}, \code{\link{extractbranch}},
  \code{\link{treeprune}} for manipulating a network using the
  branch labels.
  
  \code{\link{linnet}} for creating a network.
}
\examples{
  # make a simple tree
  m <- simplenet$m
  m[8,10] <- m[10,8] <- FALSE
  L <- linnet(vertices(simplenet), m)
  plot(L, main="")
  # compute branch labels 
  tb <- treebranchlabels(L, 1)
  tbc <- paste0("[", tb, "]")
  text(vertices(L), labels=tbc, cex=2)
}
\keyword{spatial}
\keyword{math}

