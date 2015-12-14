\name{Extract.linnet}
\alias{[.linnet}
\title{Extract Subset of Linear Network}
\description{
  Extract a subset of a linear network.
}
\usage{
  \method{[}{linnet}(x, i, \dots)
}
\arguments{
  \item{x}{
    A linear network (object of class \code{"linnet"}).
  }
  \item{i}{
    Spatial window defining the subregion.
    An object of class \code{"owin"}.
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
}
\examples{
  plot(simplenet)
  B <- owin(c(0,1),c(0.2,0.5))
  plot(simplenet[B], add=TRUE, col="red")
  plot(B, add=TRUE, border="green")
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{manip}
