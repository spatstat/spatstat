\name{rMatClust}
\alias{rMatClust}
\title{Simulate Matern Cluster Process}
\description{
  Generate a random point pattern, a simulated realisation of the
  Mat\'ern Cluster Process.
}
\usage{
 rMatClust(kappa, r, mu, win = owin(c(0,1),c(0,1)))
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{r}{
    Radius parameter of the clusters.
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
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).

  Additionally,  some intermediate results of the simulation are
  returned as attributes of this point pattern.
  See \code{\link{rNeymanScott}}.
}
\details{
  This algorithm generates a realisation of Mat\'ern's cluster process
  inside the window \code{win}. The process is constructed by first
  generating a Poisson point process of ``parent'' points 
  with intensity \code{kappa}. Then each parent point is
  replaced by a random cluster of points, the number of points in each
  cluster being random with a Poisson (\code{mu}) distribution,
  and the points being placed independently and uniformly inside
  a disc of radius \code{r} centred on the parent point.

  In this implementation, parent points are not restricted to lie in the
  window; the parent process is effectively the uniform
  Poisson process on the infinite plane.

  This classical model can be fitted to data by the method of minimum contrast,
  using \code{\link{matclust.estK}} or \code{\link{kppm}}.
  
  The algorithm can also generate spatially inhomogeneous versions of
  the Mat\'ern cluster process:
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
    with intensity function equal to the \emph{average} value of
    \code{mu} inside the disc of radius \code{r} centred on the parent
    point, and zero intensity outside this disc.
  }
  When the parents are homogeneous (\code{kappa} is a single number)
  and the offspring are inhomogeneous (\code{mu} is a
  function or pixel image), the model can be fitted to data
  using \code{\link{kppm}}, or using \code{\link{matclust.estK}}
  applied to the inhomogeneous \eqn{K} function. 
}
\seealso{
  \code{\link{rpoispp}},
  \code{\link{rThomas}},
  \code{\link{rGaussPoisson}},
  \code{\link{rNeymanScott}},
  \code{\link{matclust.estK}},
  \code{\link{kppm}}.
}
\examples{
 # homogeneous
 X <- rMatClust(10, 0.05, 4)
 # inhomogeneous
 Z <- as.im(function(x,y){ 4 * exp(2 * x - 1) }, owin())
 Y <- rMatClust(10, 0.05, Z)
}
\references{
  Mat\'ern, B. (1960)
  \emph{Spatial Variation}.
  Meddelanden fraan Statens Skogsforskningsinstitut,
  volume 59, number 5.
  Statens Skogsforskningsinstitut, Sweden.

  Mat\'ern, B. (1986)
  \emph{Spatial Variation}.
  Lecture Notes in Statistics 36, Springer-Verlag, New York.

  Waagepetersen, R. (2006)
  An estimating function approach to inference for inhomogeneous
  Neyman-Scott processes.
  Submitted for publication.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

