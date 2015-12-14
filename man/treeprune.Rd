\name{treeprune}
\alias{treeprune}
\title{
  Prune Tree to Given Level
}
\description{
  Prune a tree by removing all the branches above a given level.
}
\usage{
treeprune(X, root = 1, level = 0)
}
\arguments{
  \item{X}{
    Object of class \code{"linnet"} or \code{"lpp"}.
  }
  \item{root}{
    Index of the root vertex amongst the vertices of \code{as.linnet(X)}.
  }
  \item{level}{
    Integer specifying the level above which the tree should be pruned.
  }
}
\details{
  The object \code{X} must be either a linear network, or a derived
  object such as a point pattern on a linear network. The linear network
  must be an acyclic graph (i.e. must not contain any loops) so that it
  can be interpreted as a tree. 
  
  This function removes all vertices
  for which \code{\link{treebranchlabels}} gives a
  string more than \code{level} characters long.
}
\value{
  Object of the same kind as \code{X}.
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
  \code{\link{linnet}}
}
\examples{
  example(treebranchlabels)
  tp <- treeprune(L, 1)
  plot(tp, add=TRUE, col="blue")
}
\keyword{spatial}
\keyword{manip}

