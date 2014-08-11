\name{rpoisppx}
\alias{rpoisppx}
\title{
  Generate Poisson Point Pattern in Any Dimensions
}
\description{
  Generate a random multi-dimensional point pattern
  using the homogeneous Poisson process.
}
\usage{
rpoisppx(lambda, domain)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process.
    A single positive number.
  }
  \item{domain}{
    Multi-dimensional box in which the process should be generated.
    An object of class \code{"boxx"}.
  }
}
\details{
  This function generates a realisation
  of the homogeneous Poisson process in multi dimensions,
  with intensity \code{lambda} (points per unit volume).
  
  The realisation is generated inside the multi-dimensional region
  \code{domain} which currently must be a rectangular box (object of
  class \code{"boxx"}).
}
\value{
  The simulated multi-dimensional point pattern
  (an object of class \code{"ppx"}).
}
\note{
  The intensity \code{lambda} is the expected number of points
  \emph{per unit volume}. 
}
\seealso{
  \code{\link{runifpointx}}, 
  \code{\link{ppx}}, 
  \code{\link{boxx}}
}
\examples{
   w <- boxx(x=c(0,1), y=c(0,1), z=c(0,1), t=c(0,3))
   X <- rpoisppx(10, w)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
