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
runiflpp(n, L, nsim=1, drop=TRUE)
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
  \item{drop}{
    Logical value indicating what to do when \code{nsim=1}.
    If \code{drop=TRUE} (the default), the result is a point pattern.
    If \code{drop=FALSE}, the result is a list with one entry which is a
    point pattern.
  }
}
\details{
  This function uses \code{\link{runifpointOnLines}}
  to generate the random points.
}
\value{
  If \code{nsim = 1} and \code{drop=TRUE},
  a point pattern on the linear network,
  i.e.\ an object of class \code{"lpp"}.
  Otherwise, a list of such point patterns.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  \adrian 
}
\seealso{
 \code{\link{rlpp}} for non-uniform random points;
 \code{\link{rpoislpp}} for Poisson point process;
 
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
