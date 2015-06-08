\name{rThomas}
\alias{rThomas}
\title{Simulate Thomas Process}
\description{
  Generate a random point pattern, a realisation of the
  Thomas cluster process.
}
\usage{
  rThomas(kappa, scale, mu, win = owin(c(0,1),c(0,1)),
          nsim=1, drop=TRUE, 
          saveLambda=FALSE, expand = 4*scale, ...,
          poisthresh=1e-6)
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{scale}{
    Standard deviation of random displacement (along each coordinate axis)
    of a point from its cluster centre.
  }
  \item{mu}{
    Mean number of points per cluster (a single positive number)
    or reference intensity for the cluster points (a function or
    a pixel image).
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
  \item{saveLambda}{
    Logical. If \code{TRUE} then the random intensity corresponding to
    the simulated parent points will also be calculated and saved,
    and returns as an attribute of the point pattern.
  }
  \item{expand}{
    Numeric. Size of window expansion for generation of parent
    points. Has a sensible default.
  }
  \item{\dots}{
    Passed to \code{\link{clusterfield}} to control the image
    resolution when \code{saveLambda=TRUE} and to
    \code{\link{clusterradius}} when \code{expand} is missing.
  }
  \item{poisthresh}{
    Numerical threshold below which the model will be treated
    as a Poisson process. See Details.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"}) if \code{nsim=1},
  or a list of point patterns if \code{nsim > 1}.

  Additionally, some intermediate results of the simulation are returned
  as attributes of this point pattern (see
  \code{\link{rNeymanScott}}). Furthermore, the simulated intensity
  function is returned as an attribute \code{"Lambda"}, if
  \code{saveLambda=TRUE}.
}
\details{
  This algorithm generates a realisation of the (`modified')
  Thomas process, a special case of the Neyman-Scott process,
  inside the window \code{win}.

  In the simplest case, where \code{kappa} and \code{mu}
  are single numbers, the algorithm 
  generates a uniform Poisson point process of \dQuote{parent} points 
  with intensity \code{kappa}. Then each parent point is
  replaced by a random cluster of \dQuote{offspring} points,
  the number of points per cluster being Poisson (\code{mu})
  distributed, and their
  positions being isotropic Gaussian displacements from the
  cluster parent location. The resulting point pattern
  is a realisation of the classical
  \dQuote{stationary Thomas process} generated inside the window \code{win}.
  This point process has intensity \code{kappa * mu}.

  The algorithm can also generate spatially inhomogeneous versions of
  the Thomas process:
  \itemize{
    \item The parent points can be spatially inhomogeneous.
    If the argument \code{kappa} is a \code{function(x,y)}
    or a pixel image (object of class \code{"im"}), then it is taken
    as specifying the intensity function of an inhomogeneous Poisson
    process that generates the parent points.
    \item The offspring points can be inhomogeneous. If the
    argument \code{mu} is a \code{function(x,y)}
    or a pixel image (object of class \code{"im"}), then it is
    interpreted as the reference density for offspring points,
    in the sense of Waagepetersen (2007).
    For a given parent point, the offspring constitute a Poisson process
    with intensity function equal to \code{mu * f},
    where \code{f} is the Gaussian probability density
    centred at the parent point. Equivalently we first generate,
    for each parent point, a Poisson (\code{mumax}) random number of
    offspring (where \eqn{M} is the maximum value of \code{mu})
    with independent Gaussian displacements from the parent
    location, and then randomly thin the offspring points, with
    retention probability \code{mu/M}.
    \item Both the parent points and the offspring points can be
    spatially inhomogeneous, as described above.
  }

  Note that if \code{kappa} is a pixel image, its domain must be larger
  than the window \code{win}. This is because an offspring point inside
  \code{win} could have its parent point lying outside \code{win}.
  In order to allow this, the simulation algorithm
  first expands the original window \code{win}
  by a distance \code{expand} and generates the Poisson process of
  parent points on this larger window. If \code{kappa} is a pixel image,
  its domain must contain this larger window.

  The intensity of the Thomas process is \code{kappa * mu}
  if either \code{kappa} or \code{mu} is a single number. In the general
  case the intensity is an integral involving \code{kappa}, \code{mu}
  and \code{f}.

  The Thomas process with homogeneous parents
  (i.e. where \code{kappa} is a single number)
  can be fitted to data using \code{\link{kppm}}.
  Currently it is not possible to fit the Thomas model
  with inhomogeneous parents.

  If the pair correlation function of the model is very close
  to that of a Poisson process, deviating by less than
  \code{poisthresh}, then the model is approximately a Poisson process,
  and will be simulated as a Poisson process with intensity
  \code{kappa * mu}, using \code{\link{rpoispp}}.
  This avoids computations that would otherwise require huge amounts
  of memory.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMatClust}},
\code{\link{rCauchy}},
\code{\link{rVarGamma}},
\code{\link{rNeymanScott}},
\code{\link{rGaussPoisson}},
\code{\link{kppm}},
\code{\link{clusterfit}}.
}

\references{
  Diggle, P. J., Besag, J. and Gleaves, J. T. (1976)
  Statistical analysis of spatial point patterns by
  means of distance methods. \emph{Biometrics} \bold{32} 659--667.

  Thomas, M. (1949) A generalisation of Poisson's binomial limit for use
  in ecology. \emph{Biometrika} \bold{36}, 18--25.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\examples{
  #homogeneous
  X <- rThomas(10, 0.2, 5)
  #inhomogeneous
  Z <- as.im(function(x,y){ 5 * exp(2 * x - 1) }, owin())
  Y <- rThomas(10, 0.2, Z)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

