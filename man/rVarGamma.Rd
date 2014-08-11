\name{rVarGamma}
\alias{rVarGamma}
\title{Simulate Neyman-Scott Point Process with Variance Gamma cluster kernel}
\description{
  Generate a random point pattern, a simulated realisation of the
  Neyman-Scott process with Variance Gamma (Bessel) cluster kernel.
}
\usage{
 rVarGamma(kappa, nu.ker, omega, mu, win = owin(), 
           eps = 0.001, nu.pcf=NULL)
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{nu.ker}{
    Shape parameter for the cluster kernel. A number greater than -1.
  }
  \item{omega}{
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
  \item{eps}{
    Threshold below which the values of the cluster kernel
    will be treated as zero for simulation purposes. 
  }
  \item{nu.pcf}{
    Alternative specifier of the shape parameter. See Details.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).

  Additionally,  some intermediate results of the simulation are
  returned as attributes of this point pattern.
  See \code{\link{rNeymanScott}}.
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
  index \code{nu.ker}. This is the parameter
  \eqn{\nu^\prime = \alpha/2-1}{nu' = alpha/2 - 1} appearing in
  equation (12) on page 126 of Jalilian et al (2013).
  Instead of specifying \code{nu.ker} the user can specify
  \code{nu.pcf} which is the parameter \eqn{\nu=\alpha-1}{nu = alpha-1}
  appearing in equation (13), page 127 of Jalilian et al (2013).
  These are related by \code{nu.pcf = 2 * nu.ker + 1}
  and \code{nu.ker = (nu.pcf - 1)/2}.
  Exactly one of \code{nu.ker} or \code{nu.pcf} must be specified.

  The scale of the kernel is determined by the argument \code{omega},
  which is the parameter
  \eqn{\eta}{eta} appearing in equations (12) and (13) of
  Jalilian et al (2013). 
  It is expressed in units of length (the same as the unit of length for 
  the window \code{win}).
  
  In this implementation, parent points are not restricted to lie in the
  window; the parent process is effectively the uniform
  Poisson process on the infinite plane.

  This model can be fitted to data by the method of minimum contrast,
  using \code{\link{cauchy.estK}}, \code{\link{cauchy.estpcf}}
  or \code{\link{kppm}}. It can also be fitted by maximum composite
  likelihood using \code{\link{kppm}}.
  
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
  using \code{\link{kppm}}, or using \code{\link{cauchy.estK}}
  or \code{\link{cauchy.estpcf}}
  applied to the inhomogeneous \eqn{K} function. 
}
\seealso{
  \code{\link{rpoispp}},
  \code{\link{rNeymanScott}},
  \code{\link{cauchy.estK}},
  \code{\link{cauchy.estpcf}},
  \code{\link{kppm}}.
}
\examples{
 # homogeneous
 X <- rVarGamma(30, 2, 0.02, 5)
 # inhomogeneous
 Z <- as.im(function(x,y){ exp(2 - 3 * x) }, W= owin())
 Y <- rVarGamma(30, 2, 0.02, Z)
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

