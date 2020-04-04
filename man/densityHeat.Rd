\name{densityHeat}
\alias{densityHeat}
\title{
  Kernel Density on a Network using Heat Equation
}
\description{
  Computes a kernel density estimate on a linear network
  using the heat equation.
}
\usage{
  densityHeat(x, sigma, \dots, weights = NULL,
              dx = NULL, dt = NULL, iterMax = 1e+06, verbose=FALSE)
}
\arguments{
  \item{x}{
    Point pattern on a linear network (object of class \code{"lpp"})
    to be smoothed.
  }
  \item{sigma}{
    Smoothing bandwidth (standard deviation of the kernel)
    in the same units as the spatial coordinates of \code{x}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} determining the
    resolution of the result. (Any other arguments are ignored.)
  }
  \item{weights}{
    Optional. Numeric vector of weights associated with the
    points of \code{x}. Weights may be positive, negative or zero.
  }
  \item{dx}{
    Optional. Spacing of the sampling points along the network.
    A single number giving a distance value in the same units as \code{x}.
  }
  \item{dt}{
    Optional. Time step in the heat equation solver.
    A single number. 
  }
  \item{iterMax}{
    Maximum number of iterations.
  }
  \item{verbose}{
    Logical value specifying whether to print progress reports.
  }
}
\details{
  Kernel smoothing is applied to the points of \code{x}
  using a kernel based on path distances in the network.
  The result is a pixel image on the linear network (class
  \code{"linim"}) which can be plotted.

  Smoothing is performed using the \dQuote{equal-split continuous} rule described in
  Section 9.2.3 of Okabe and Sugihara (2012).
  However, the actual computation is performed rapidly, by solving the classical
  time-dependent heat equation on the network, as described in McSwiggan et al (2016).
  Computational time is short, but increases quadratically with \code{sigma}.
}
\section{Infinite bandwidth}{
  If \code{sigma=Inf}, the resulting density estimate is 
  constant over all locations,
  and is equal to the average density of points per unit length.
  (If the network is not connected, then this rule
  is applied separately to each connected component of the network).
}
\value{
  A pixel image on the linear network (object of class \code{"linim"}).
}
\references{
  McSwiggan, G., Baddeley, A. and Nair, G. (2016)
  Kernel density estimation on a linear network.
  \emph{Scandinavian Journal of Statistics} \bold{44}, 324--345.

  Okabe, A. and Sugihara, K. (2012)
  \emph{Spatial analysis along networks}.
  Wiley.
}
\author{
  \adrian and Greg McSwiggan.
}
\seealso{
  \code{\link{density.lpp}}
}
\examples{
  X <- runiflpp(3, simplenet)
  D <- densityHeat(X, 0.2, verbose=FALSE)
  plot(D, style="w", main="", adjust=2)
  Dw <- densityHeat(X, 0.2, weights=c(1,2,-1), verbose=FALSE)
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
