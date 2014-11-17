\name{runiflpp}
\alias{runiflpp}
\title{
  Uniform Random Points on a Linear Network
}
\description{
  Generates \eqn{n} random points, independently and
  uniformly distributed, on a linear network.
}
\usage{
runiflpp(n, L, nsim=1)
}
\arguments{
  \item{n}{
    Number of random points to generate.
    A nonnegative integer, or a vector of integers
    specifying the number of points of each type.
  }
  \item{L}{
    A linear network (object of class \code{"linnet"},
    see \code{\link{linnet}}).
  }
  \item{nsim}{Number of simulated realisations to generate.}
}
\details{
  This function uses \code{\link{runifpointOnLines}}
  to generate the random points.
}
\value{
  If \code{nsim = 1}, a point pattern on the linear network,
  i.e.\ an object of class \code{"lpp"}.
  If \code{nsim > 1}, a list of such point patterns.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
 \code{\link{rpoislpp}},
 \code{\link{lpp}},
 \code{\link{linnet}}
}
\examples{
   data(simplenet)
   X <- runiflpp(10, simplenet)
   plot(X)
  # marked
   Z <- runiflpp(c(a=10, b=3), simplenet)
}
\keyword{spatial}
\keyword{datagen}
