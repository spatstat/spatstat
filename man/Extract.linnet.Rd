\name{Extract.linnet}
\alias{[.linnet}
\title{Extract Subset of Linear Network}
\description{
  Extract a subset of a linear network.
}
\usage{
  \method{[}{linnet}(x, i, \dots, snip=TRUE)
}
\arguments{
  \item{x}{
    A linear network (object of class \code{"linnet"}).
  }
  \item{i}{
    Spatial window defining the subregion.
    An object of class \code{"owin"}.
  }
  \item{snip}{
    Logical. If \code{TRUE} (the default), segments of \code{x}
    which cross the boundary of \code{i} will be cut by the boundary.
    If \code{FALSE}, these segments will be deleted.
  }
  \item{\dots}{Ignored.}
}
\value{
  Another linear network (object of class \code{"linnet"}).
}
\details{
  This function computes the intersection
  between the linear network \code{x} and the domain specified by \code{i}.
  
  This function is a method for the subset operator \code{"["} for
  linear networks (objects of class \code{"linnet"}). It is provided
  mainly for completeness.

  The index \code{i} should be a window.

  The argument \code{snip} specifies what to do with segments of \code{x}
  which cross the boundary of \code{i}. 
  If \code{snip=FALSE}, such segments are simply deleted.  
  If \code{snip=TRUE} (the default), such segments are cut into pieces by the
  boundary of \code{i}, and those pieces which lie inside the window
  \code{i} are included in the resulting network. 
}
\examples{
  p <- par(mfrow=c(1,2), mar=0.2+c(0,0,1,0))
  B <- owin(c(0.1,0.7),c(0.19,0.6))

  plot(simplenet, main="x[w, snip=TRUE]")
  plot(simplenet[B], add=TRUE, col="green", lwd=3)
  plot(B, add=TRUE, border="red", lty=3)

  plot(simplenet, main="x[w, snip=FALSE]")
  plot(simplenet[B, snip=FALSE], add=TRUE, col="green", lwd=3)
  plot(B, add=TRUE, border="red", lty=3)

  par(p)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
  Rolf Turner \email{r.turner@auckland.ac.nz},
  Ege Rubak \email{rubak@math.aau.dk}
  and Suman Rakshit.
}
\keyword{spatial}
\keyword{manip}
