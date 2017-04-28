\name{as.linfun} 
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

  \method{as.linfun}{lintess}(X, \dots, values, navalue=NA)
}
\arguments{
  \item{X}{
    Some kind of data to be converted.
  }
  \item{\dots}{
    Other arguments passed to methods.
  }
  \item{values}{
    Optional. Vector of function values,
    one entry associated with each tile of the tessellation.
  }
  \item{navalue}{
    Optional. Function value associated with locations that
    do not belong to a tile of the tessellation.
  }
}
\details{
  An object of class \code{"linfun"} represents
  a function defined on a linear network.

  The function \code{as.linfun} is generic. The method \code{as.linfun.linim}
  converts objects of class \code{"linim"} (pixel images on a linear
  network) to functions on the network.

  The method \code{as.linfun.lintess} converts a tessellation
  on a linear network into a function with a different value
  on each tile of the tessellation.
  If the argument \code{values} is missing or null,
  then the function returns factor values identifying which tile
  contains each given point. If \code{values} is given,
  it should be a vector with one entry for each tile of the tessellation:
  any point lying in tile number \code{i} will return the value \code{v[i]}.
}
\value{
  Object of class \code{"linfun"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{linfun}}
}
\examples{
   X <- runiflpp(2, simplenet)
   Y <- runiflpp(5, simplenet)

   # image on network
   D <- density(Y, 0.1, verbose=FALSE)

   f <- as.linfun(D)
   f
   f(X)

   # tessellation on network
   Z <- lineardirichlet(Y)
   g <- as.linfun(Z)
   g(X)
   h <- as.linfun(Z, values = runif(5))
   h(X)
}
\keyword{spatial}
\keyword{manip}
