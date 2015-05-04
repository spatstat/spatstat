\name{rpoispp3}
\alias{rpoispp3}
\title{
  Generate Poisson Point Pattern in Three Dimensions
}
\description{
  Generate a random three-dimensional point pattern
  using the homogeneous Poisson process.
}
\usage{
rpoispp3(lambda, domain = box3(), nsim=1, drop=TRUE)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process.
    A single positive number.
  }
  \item{domain}{
    Three-dimensional box in which the process should be generated.
    An object of class \code{"box3"}.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  If \code{nsim = 1} and \code{drop=TRUE}, a point pattern in
  three dimensions (an object of class \code{"pp3"}).
  If \code{nsim > 1}, a list of such point patterns.
}
\details{
  This function generates a realisation
  of the homogeneous Poisson process in three dimensions,
  with intensity \code{lambda} (points per unit volume).
  
  The realisation is generated inside the three-dimensional region
  \code{domain} which currently must be a rectangular box (object of
  class \code{"box3"}).
}
\note{
  The intensity \code{lambda} is the expected number of points
  \emph{per unit volume}. 
}
\seealso{
  \code{\link{runifpoint3}}, 
  \code{\link{pp3}}, 
  \code{\link{box3}}
}
\examples{
   X <- rpoispp3(50)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
