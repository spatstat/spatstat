\name{as.linfun} %DontDeclareMethods
\alias{as.linfun}
\alias{as.linfun.linim}
\alias{as.linfun.lintess}
\title{
  Convert Data to a Function on a Linear Network
}
\description{
  Convert some kind of data to an object of class \code{"linfun"}
  representing a function on a linear network.
}
\usage{
  as.linfun(X, \dots)

  \method{as.linfun}{linim}(X, \dots)

  \method{as.linfun}{lintess}(X, \dots)
}
\arguments{
  \item{X}{
    Some kind of data to be converted.
  }
  \item{\dots}{
    Other arguments passed to methods.
  }
}
\details{
  An object of class \code{"linfun"} represents
  a function defined on a linear network.

  The function \code{as.linfun} is generic. The method \code{as.linfun.linim}
  converts objects of class \code{"linim"} (pixel images on a linear
  network) to functions on the network.
  The method \code{as.linfun.lintess} converts a tessellation
  on a linear network into a function identifying which tile contains
  any given point.
}
\value{
  Object of class \code{"linfun"}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\seealso{
  \code{\link{linfun}}
}
\examples{
   D <- density(runiflpp(4, simplenet), 0.1, verbose=FALSE)

   f <- as.linfun(D)
   f
   
   X <- runiflpp(2, simplenet)
   f(X)
}
\keyword{spatial}
\keyword{manip}
