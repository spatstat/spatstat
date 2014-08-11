\name{rThomas}
\alias{rThomas}
\title{Simulate Thomas Process}
\description{
  Generate a random point pattern, a realisation of the
  Thomas cluster process.
}
\usage{
 rThomas(kappa, sigma, mu, win = owin(c(0,1),c(0,1)))
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number.
  }
  \item{sigma}{
    Standard deviation of random displacement (along each coordinate axis)
    of a point from its cluster centre.
  }
  \item{mu}{
    Expected number of points per cluster.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).

  Additionally,  some intermediate results of the simulation are
  returned as attributes of this point pattern.
  See \code{\link{rNeymanScott}}.
}
\details{
  This algorithm generates a realisation of the
  Thomas process, a special case of the Neyman-Scott process.
  
  The algorithm 
  generates a uniform Poisson point process of ``parent'' points 
  with intensity \code{kappa}. Then each parent point is
  replaced by a random cluster of points, the number of points
  per cluster being Poisson (\code{mu}) distributed, and their
  positions being isotropic Gaussian displacements from the
  cluster parent location.

  This classical model can be fitted to data by the method of minimum contrast,
  using \code{\link{thomas.estK}} or \code{\link{kppm}}.
  
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
    in the sense of Waagepetersen (2006).
    For a given parent point, the offspring constitute a Poisson process
    with intensity function equal to \code{mu(x,y) * f(x,y)}
    where \code{f} is the Gaussian density centred at the parent point.
  }
  When the parents are homogeneous (\code{kappa} is a single number)
  and the offspring are inhomogeneous (\code{mu} is a
  function or pixel image), the model can be fitted to data
  using \code{\link{kppm}}, or 
  using \code{\link{thomas.estK}} applied to the inhomogeneous
  \eqn{K} function.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMatClust}},
\code{\link{rGaussPoisson}},
\code{\link{rNeymanScott}},
\code{\link{thomas.estK}},
\code{\link{kppm}}
}
\references{
  Waagepetersen, R. (2006)
  An estimating function approach to inference for inhomogeneous
  Neyman-Scott processes.
  Submitted for publication.
}
\examples{
  #homogeneous
  X <- rThomas(10, 0.2, 5)
  #inhomogeneous
  Z <- as.im(function(x,y){ 5 * exp(2 * x - 1) }, owin())
  Y <- rThomas(10, 0.2, Z)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

