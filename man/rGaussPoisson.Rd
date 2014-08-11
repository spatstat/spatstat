\name{rGaussPoisson}
\alias{rGaussPoisson}
\title{Simulate Gauss-Poisson Process}
\description{
  Generate a random point pattern, a simulated realisation of the
  Gauss-Poisson Process.
}
\usage{
 rGaussPoisson(kappa, r, p2, win = owin(c(0,1),c(0,1)))
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{r}{
    Diameter of each cluster that consists of exactly 2 points.
  }
  \item{p2}{
    Probability that a cluster contains exactly 2 points.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).

  Additionally, some intermediate results of the simulation are
  returned as attributes of this point pattern.
  See \code{\link{rNeymanScott}}.
}
\details{
  This algorithm generates a realisation of the Gauss-Poisson
  point process inside the window \code{win}.
  The process is constructed by first
  generating a Poisson point process of parent points 
  with intensity \code{kappa}. Then each parent point is either retained
  (with probability \code{1 - p2})
  or replaced by a pair of points at a fixed distance \code{r} apart
  (with probability \code{p2}). In the case of clusters of 2 points,
  the line joining the two points has uniform random orientation.

  In this implementation, parent points are not restricted to lie in the
  window; the parent process is effectively the uniform
  Poisson process on the infinite plane.
}
\seealso{
  \code{\link{rpoispp}},
  \code{\link{rThomas}},
  \code{\link{rMatClust}},
  \code{\link{rNeymanScott}}
}
\examples{
 pp <- rGaussPoisson(30, 0.07, 0.5)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

