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
  clusterfield(model, locations = NULL, \dots)

  \method{clusterfield}{character}(model, locations = NULL, \dots)

  \method{clusterfield}{function}(model, locations = NULL, \dots, mu = NULL)

  \method{clusterfield}{kppm}(model, locations = NULL, \dots)
}
\arguments{
  \item{model}{
    Cluster model. Either a fitted cluster model (object of class
    \code{"kppm"}), a character string specifying the type of cluster
    model, or a function defining the cluster kernel. See Details.
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
  \item{mu}{
     Mean number of offspring per cluster. 
     A single number or a pixel image.
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

  The function method also accepts the optional parameter \code{mu}
  (defaulting to 1) specifying the mean number of points per cluster (as
  a numeric) or the inhomogeneous reference cluster intensity (as an
  \code{"im"} object or a \code{function(x,y)}). The interpretation of
  \code{mu} is as explained in the simulation functions referenced in
  the See Also section below.

  For the character method \code{model} must be one of:
  \code{model="Thomas"} for the Thomas process,
  \code{model="MatClust"} for the Matern cluster process,
  \code{model="Cauchy"} for the Neyman-Scott cluster process with
  Cauchy kernel, or \code{model="VarGamma"} for the Neyman-Scott
  cluster process with Variance Gamma kernel. For all these models the
  parameter \code{scale} is required and passed through \code{\dots} as
  well as the parameter \code{nu} when \code{model="VarGamma"}. This
  method calls \code{clusterfield.function} so the parameter \code{mu}
  may also be passed through \code{\dots} and will be interpreted as
  explained above.

  The kppm method extracts the relevant information from the fitted
  model (including \code{mu}) and calls \code{clusterfield.function}.
  
}
\value{
  A pixel image (object of class \code{"im"}).
}
\seealso{
  \code{\link{density.ppp}} and \code{\link{kppm}}

  Simulation algorithms for cluster models:
  \code{\link{rCauchy}}
  \code{\link{rMatClust}}
  \code{\link{rThomas}}
  \code{\link{rVarGamma}}
}
\examples{
  # method for fitted model
  fit <- kppm(redwood~1, "Thomas")
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
