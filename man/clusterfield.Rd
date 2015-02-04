\name{clusterfield}
\alias{clusterfield}
\alias{clusterfield.character}
\alias{clusterfield.function}
\alias{clusterfield.kppm}
\title{Field of clusters}
\description{
  Calculate the superposition of cluster kernels at the location of a
  point pattern.
}
\usage{
  clusterfield(clusters, locations = NULL, \dots)

  \method{clusterfield}{function}(clusters, locations = NULL, \dots)

  \method{clusterfield}{kppm}(clusters, locations = NULL, \dots)
}
\arguments{
  \item{clusters}{
    A function of defining the cluster kernel or a fitted cluster model
    (object of class \code{"kppm"}) from which it can be extracted or a
    string specifying a model. See Details.
  }
  \item{locations}{
    A point pattern giving the locations of the kernels. Defaults to the
    centroid of the observation window for the \code{"kppm"} method and
    to the center of a unit square otherwise.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{density.ppp}} or the
    cluster kernel. See Details.
  }
}
\details{
  The actual calculations are preformed by \code{\link{density.ppp}} and
  \code{\dots} arguments are passed thereto for control over the pixel
  resolution etc. (These arguments are then passed on to \code{\link{pixellate.ppp}}
  and \code{\link{as.mask}}.)

  For the function method the given kernel function should accept
  vectors of x and y coordinates as its first two arguments. Any
  additional arguments may be passed through the \code{\dots}.

  For the character method \code{clusters} must be one of:
  \code{clusters="Thomas"} for the Thomas process,
  \code{clusters="MatClust"} for the Matern cluster process,
  \code{clusters="Cauchy"} for the Neyman-Scott cluster process
  with Cauchy kernel, or
  \code{clusters="VarGamma"} for the Neyman-Scott cluster process
  with Variance Gamma kernel.
  For all these models the parameter \code{scale} is required and passed
  through the dots as well as the parameter \code{nu} when
  \code{clusters="VarGamma"}.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\seealso{
  \code{\link{density.ppp}} and \code{\link{kppm}}
}
\examples{
  # method for fitted model
  fit <- kppm(redwood, ~1, "Thomas")
  clusterfield(fit, eps = 0.01)

  # method for functions
  kernel <- function(x,y,scal) { 
      r <- sqrt(x^2 + y^2)
      ifelse(r > 0,
             dgamma(r, shape=5, scale=scal)/(2 * pi * r),
             0)               
  }
  X <- runifpoint(10)
  clusterfield(kernel, X, scal=0.05)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}.
}
\keyword{spatial}
