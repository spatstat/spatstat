\name{rVarGamma}
\alias{rVarGamma}
\title{Simulate Neyman-Scott Point Process with Variance Gamma cluster kernel}
\description{
  Generate a random point pattern, a simulated realisation of the
  Neyman-Scott process with Variance Gamma (Bessel) cluster kernel.
}
\usage{
 rVarGamma(kappa, nu, scale, mu, win = owin(), 
           thresh = 0.001, nsim=1, drop=TRUE,
           saveLambda=FALSE, expand = NULL, ...,
           poisthresh=1e-6)
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{nu}{
    Shape parameter for the cluster kernel. A number greater than -1.
  }
  \item{scale}{
    Scale parameter for cluster kernel. Determines the size of clusters.
    A positive number in the same units as the spatial coordinates.
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
  \item{thresh}{
    Threshold relative to the cluster kernel value at the origin (parent
    location) determining when the cluster kernel will be treated as
    zero for simulation purposes. Will be overridden by argument
    \code{expand} if that is given.
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
    points. By default determined by calling
    \code{\link{clusterradius}} with the numeric threshold value given
    in \code{thresh}.
  }
  \item{\dots}{
    Passed to \code{\link{clusterfield}} to control the image resolution
    when \code{saveLambda=TRUE} and to \code{\link{clusterradius}} when
    \code{expand} is missing or \code{NULL}.
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
  This algorithm generates a realisation of the Neyman-Scott process
  with Variance Gamma (Bessel) cluster kernel, inside the window \code{win}.

  The process is constructed by first
  generating a Poisson point process of ``parent'' points 
  with intensity \code{kappa}. Then each parent point is
  replaced by a random cluster of points, the number of points in each
  cluster being random with a Poisson (\code{mu}) distribution,
  and the points being placed independently and uniformly
  according to a Variance Gamma kernel.

  The shape of the kernel is determined by the dimensionless
  index \code{nu}. This is the parameter
  \eqn{\nu^\prime = \alpha/2-1}{nu' = alpha/2 - 1} appearing in
  equation (12) on page 126 of Jalilian et al (2013).

  The scale of the kernel is determined by the argument \code{scale},
  which is the parameter
  \eqn{\eta}{eta} appearing in equations (12) and (13) of
  Jalilian et al (2013). 
  It is expressed in units of length (the same as the unit of length for 
  the window \code{win}).
  
  In this implementation, parent points are not restricted to lie in the
  window; the parent process is effectively the uniform
  Poisson process on the infinite plane.

  This model can be fitted to data by the method of minimum contrast,
  maximum composite likelihood or Palm likelihood using
  \code{\link{kppm}}.
  
  The algorithm can also generate spatially inhomogeneous versions of
  the cluster process:
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
    in the sense of Waagepetersen (2006).
  }
  When the parents are homogeneous (\code{kappa} is a single number)
  and the offspring are inhomogeneous (\code{mu} is a
  function or pixel image), the model can be fitted to data
  using \code{\link{kppm}}, or using \code{\link{vargamma.estK}}
  or \code{\link{vargamma.estpcf}}
  applied to the inhomogeneous \eqn{K} function. 

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
  \code{\link{rNeymanScott}},
  \code{\link{kppm}}.

  \code{\link{vargamma.estK}},
  \code{\link{vargamma.estpcf}}.
}
\examples{
 # homogeneous
 X <- rVarGamma(30, 2, 0.02, 5)
 # inhomogeneous
 ff <- function(x,y){ exp(2 - 3 * abs(x)) }
 Z <- as.im(ff, W= owin())
 Y <- rVarGamma(30, 2, 0.02, Z)
 YY <- rVarGamma(ff, 2, 0.02, 3)
}
\references{
  Jalilian, A., Guan, Y. and Waagepetersen, R. (2013)
  Decomposition of variance for spatial Cox processes.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 119-137.
  
  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{Abdollah Jalilian and Rasmus Waagepetersen.
  Adapted for \pkg{spatstat} by Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}
\keyword{datagen}

