\name{connected.linnet}
\alias{connected.linnet}
\title{
  Connected Components of a Linear Network
}
\description{
  Find the topologically-connected components of a linear network.
}
\usage{
\method{connected}{linnet}(X, \dots, what = c("labels", "components"))
}
\arguments{
  \item{X}{
    A linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{what}{
    Character string specifying the kind of result.
  }
}
\details{
  The function \code{connected} is generic. This is the method for
  linear networks (objects of class \code{"linnet"}).

  Two vertices of the network are connected if they are joined by a path
  in the network. This function divides the network into subsets, such
  that all points in a subset are connected to each other.

  If \code{what="labels"} the return value is a factor with one entry
  for each vertex of \code{X}, identifying which connected component the
  vertex belongs to. 

  If \code{what="components"} the return value is a list of linear
  networks, which are the connected components of \code{X}.
}
\value{
  If \code{what="labels"}, a factor.
  If \code{what="components"}, a list of linear networks.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  and Suman Rakshit.
}
\seealso{
  \code{\link{thinNetwork}}
}
\examples{
   # remove some edges from a network to make it disconnected
   plot(simplenet, col="grey", main="", lty=2)
   A <- thinNetwork(simplenet, retainedges=-c(3,5))
   plot(A, add=TRUE, lwd=2)
   # find the connected components
   connected(A)
   cA <- connected(A, what="components")
   plot(cA[[1]], add=TRUE, col="green", lwd=2)
   plot(cA[[2]], add=TRUE, col="blue", lwd=2)
}
\keyword{spatial}
\keyword{manip}
