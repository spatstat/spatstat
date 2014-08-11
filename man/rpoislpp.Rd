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
rpoislpp(lambda, L, ...)
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
  }
  \item{\dots}{
    Arguments passed to \code{\link{rpoisppOnLines}}.
  }
}
\details{
  This function uses \code{\link{rpoisppOnLines}}
  to generate the random points.
}
\value{
  A point pattern on the linear network,
  i.e.\ an object of class \code{"lpp"}.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
 \code{\link{runiflpp}},
 \code{\link{lpp}},
 \code{\link{linnet}}
}
\examples{
   data(simplenet)
   X <- rpoislpp(5, simplenet)
   plot(X)
   # multitype
   X <- rpoislpp(c(a=5, b=5), simplenet)
}
\keyword{spatial}
\keyword{datagen}
