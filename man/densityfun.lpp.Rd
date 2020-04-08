\name{densityfun.lpp}
\alias{densityfun.lpp}
\title{
  Kernel Estimate of Intensity on a Linear Network as a Spatial Function
}
\description{
  Computes a kernel estimate of the intensity of a point process on a linear
  network, and returns the intensity estimate as a function of spatial location.
}
\usage{
\method{densityfun}{lpp}(X, sigma, \dots, weights=NULL, nsigma=1, verbose=FALSE)
}
\arguments{
  \item{X}{
    Point pattern on a linear network
    (object of class \code{"lpp"}).
  }
  \item{sigma}{
    Bandwidth of kernel (standard deviation of Gaussian kernel),
    in the same units of length as \code{X}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.lpp}} to control
    the discretisation.
  }
  \item{weights}{
    Optional numeric vector of weights associated
    with the points of \code{X}.
  }
  \item{nsigma}{
    Integer. The number of different bandwidths for which a result
    should be returned.
    If \code{nsigma=1} (the default), the result is a function giving
    kernel estimate with bandwidth \code{sigma}. If \code{nsigma > 1},
    the result is a function with an additional argument \code{k}
    containing the kernel estimates for the \code{nsigma+1}
    equally-spaced time steps from \code{0} to \code{sigma^2}.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  Kernel smoothing is applied to the points of \code{X}
  using the diffusion algorithm of McSwiggan et al (2016).
  The result is a function on the linear network
  (object of class \code{"linfun"}) that can be printed, plotted
  and evaluated at any location. 
  
  This is a method for the generic function
  \code{\link[spatstat]{densityfun}} for the class
  \code{"lpp"} of point patterns on a linear network.
}
\value{
  Function on a linear network (object of class \code{"linfun"}).

  If \code{nsigma=1} (the default), the result is a function giving
  kernel estimate with bandwidth \code{sigma}.

  If \code{nsigma > 1}, the result is a function
  with an additional argument \code{k}. If \code{k} is specified,
  the function returns the kernel estimate for
  bandwidth \code{tau = sigma * sqrt(k/nsigma)}.
  If \code{k} is not specified, results are returned
  for all \code{k = 1, 2, ..., nsigma}.
  
  The result also has attributes
  \itemize{
    \item \code{attr(result, "dt")} giving the
    time step \eqn{\Delta t}{Delta t};
    \item \code{attr(result, "dx")} giving the spacing \eqn{\Delta x}{Delta x}
    between sample points in the numerical algorithm;
    \item \code{attr(result, "sigma")}
    giving the smoothing bandwidth \eqn{\sigma}{sigma} used
    (or the successive bandwidths used at each sampled time step,
    if \code{nsigma > 1}).
  }
}
\references{
  McSwiggan, G., Baddeley, A. and Nair, G. (2016)
  Kernel Density Estimation on a Linear Network.
  \emph{Scandinavian Journal of Statistics} \bold{44}, 324--345.
}
\author{
  Greg McSwiggan, with tweaks by \adrian.
}
\seealso{
  \code{\link[spatstat]{density.lpp}} which returns a pixel image
  on the linear network.
  
  \code{\link[spatstat]{methods.linfun}} for methods applicable to
  \code{"linfun"} objects.
}
\examples{
  X <- unmark(chicago)
  # single bandwidth
  g <- densityfun(X, 30)
  plot(g)
  Y <- X[1:5]
  g(Y)
  # weighted
  gw <- densityfun(X, 30, weights=runif(npoints(X)))
  # sequence of bandwidths 
  g10 <- densityfun(X, 30, nsigma=10)
  g10(Y, k=10)
  g10(Y)
  plot(as.linim(g10, k=5))
}
\keyword{spatial}
\keyword{nonparametric}
