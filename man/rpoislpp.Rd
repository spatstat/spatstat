\name{rpoislpp}
\alias{rpoislpp}
\title{
  Poisson Point Process on a Linear Network
}
\description{
  Generates a realisation of the Poisson point process
  with specified intensity on the given linear network.
}
\usage{
rpoislpp(lambda, L, ..., nsim=1)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process. 
    A single number, a \code{function(x,y)}, a pixel image
    (object of class \code{"im"}), or a vector of numbers,
    a list of functions, or a list of images.
  }
  \item{L}{
    A linear network (object of class \code{"linnet"},
    see \code{\link{linnet}}).
    Can be omitted in some cases: see Details.
  }
  \item{\dots}{
    Arguments passed to \code{\link{rpoisppOnLines}}.
  }
  \item{nsim}{Number of simulated realisations to generate.}
}
\details{
  This function uses \code{\link{rpoisppOnLines}}
  to generate the random points.

  Argument \code{L} can be omitted, and defaults to \code{as.linnet(lambda)},
  when \code{lambda} is a function on a linear network (class
  \code{"linfun"}) or a pixel image on a linear network (\code{"linim"}).
}
\value{
  If \code{nsim = 1}, a point pattern on the linear network,
  i.e.\ an object of class \code{"lpp"}.
  If \code{nsim > 1}, a list of such point patterns.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  \adrian 
  
}
\seealso{
 \code{\link{runiflpp}},
 \code{\link{lpp}},
 \code{\link{linnet}}
}
\examples{
   X <- rpoislpp(5, simplenet)
   plot(X)
   # multitype
   X <- rpoislpp(c(a=5, b=5), simplenet)
}
\keyword{spatial}
\keyword{datagen}
