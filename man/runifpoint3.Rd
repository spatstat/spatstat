\name{runifpoint3}
\alias{runifpoint3}
\title{
  Generate N Uniform Random Points in Three Dimensions
}
\description{
  Generate a random point pattern containing \code{n}
  independent, uniform random points in three dimensions.
}
\usage{
runifpoint3(n, domain = box3(), nsim=1)
}
\arguments{
  \item{n}{
    Number of points to be generated.
  }
  \item{domain}{
    Three-dimensional box in which the process should be generated.
    An object of class \code{"box3"}.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
}
\value{
  If \code{nsim = 1}, a point pattern in
  three dimensions (an object of class \code{"pp3"}).
  If \code{nsim > 1}, a list of such point patterns.
}
\details{
  This function generates \code{n} independent random points,
  uniformly distributed in the three-dimensional box \code{domain}.
}
\seealso{
  \code{\link{rpoispp3}}, 
  \code{\link{pp3}}, 
  \code{\link{box3}}
}
\examples{
   X <- runifpoint3(50)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
