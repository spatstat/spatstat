\name{deletebranch}
\alias{deletebranch}
\alias{deletebranch.linnet}
\alias{deletebranch.lpp}
\alias{extractbranch}
\alias{extractbranch.linnet}
\alias{extractbranch.lpp}
\title{
  Delete or Extract a Branch of a Tree
}
\description{
  Deletes or extracts a given branch of a tree.
}
\usage{
deletebranch(X, \dots)

\method{deletebranch}{linnet}(X, code, labels, \dots)

\method{deletebranch}{lpp}(X, code, labels, \dots)

extractbranch(X, \dots)

\method{extractbranch}{linnet}(X, code, labels, \dots, which=NULL)

\method{extractbranch}{lpp}(X, code, labels, \dots, which=NULL)
}
\arguments{
  \item{X}{
    Linear network (object of class \code{"linnet"})
    or point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{code}{
    Character string. Label of the branch to be deleted or extracted.
  }
  \item{labels}{
    Vector of character strings. Branch labels for the
    vertices of the network, usually obtained from
    \code{\link{treebranchlabels}}.
  }
  \item{\dots}{Arguments passed to methods.}
  \item{which}{
    Logical vector indicating which vertices of the network
    should be extracted. Overrides \code{code} and \code{labels}.
  }
}
\details{
  The linear network \code{L <- X} or \code{L <- as.linnet(X)}
  must be a tree, that is, it has no loops.

  The argument \code{labels} should be a character vector
  giving tree branch labels for each vertex of the network.
  It is usually obtained by calling \code{\link{treebranchlabels}}.

  The branch designated by the string \code{code} will be deleted
  or extracted.

  The return value is the result of deleting or extracting
  this branch from \code{X}
  along with any data associated with this branch (such as points or marks).
}
\value{
  Another object of the same type as \code{X}
  obtained by deleting or extracting the specified branch.
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
\code{\link{treebranchlabels}},
\code{\link{branchlabelfun}},
\code{\link{linnet}}
}
\examples{
  example(treebranchlabels)
  LminusB <- deletebranch(L, "b", tb)
  plot(LminusB, add=TRUE, col="green")
  LB <- extractbranch(L, "b", tb)
  plot(LB, add=TRUE, col="red")
}
\keyword{spatial}
\keyword{manip}

