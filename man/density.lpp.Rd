\name{density.lpp}
\alias{density.lpp}
\alias{density.splitppx}
\title{
  Kernel Estimate of Intensity on a Linear Network
}
\description{
  Estimates the intensity of a point process on a linear network
  by applying kernel smoothing to the point pattern data.
}
\usage{
\method{density}{lpp}(x, sigma=NULL, \dots,
        weights=NULL,
        distance=c("path", "euclidean"),
        continuous=TRUE,
        kernel="gaussian")

\method{density}{splitppx}(x, sigma=NULL, \dots)
}
\arguments{
  \item{x}{
    Point pattern on a linear network (object of class \code{"lpp"})
    to be smoothed.
  }
  \item{sigma}{
    Smoothing bandwidth (standard deviation of the kernel).
    A single numerical value
    in the same units as the spatial coordinates of \code{x}.
  }
  \item{\dots}{
    Additional arguments controlling the algorithm
    and the spatial resolution of the result.
    These arguments are passed either to
    \code{\link{densityQuick.lpp}}, \code{\link{densityHeat}}
    or \code{\link{densityEqualSplit}} depending on the algorithm chosen.
  }
  \item{weights}{
    Optional. Numeric vector of weights associated with the
    points of \code{x}. Weights may be positive, negative or zero.
  }
  \item{distance}{
    Character string (partially matched) specifying whether to use
    a kernel based on paths in the network (\code{distance="path"}, the default)
    or a two-dimensional kernel (\code{distance="euclidean"}). 
  }
  \item{kernel}{
    Character string specifying the smoothing kernel.
    See \code{\link{dkernel}} for possible options.
  }
  \item{continuous}{
    Logical value indicating whether to compute the
    \dQuote{equal-split continuous} smoother (\code{continuous=TRUE}, the
    default) or the \dQuote{equal-split discontinuous} smoother
    (\code{continuous=FALSE}). Applies only when \code{distance="path"}.
  }
}
\details{
  Kernel smoothing is applied to the points of \code{x}
  using either a kernel based on path distances in the network,
  or a two-dimensional kernel.
  The result is a pixel image on the linear network (class
  \code{"linim"}) which can be plotted.

  \itemize{
    \item If \code{distance="path"} (the default)
    then the smoothing is performed
    using a kernel based on path distances in the network, as described in
    described in Okabe and Sugihara (2012) and McSwiggan et al (2016).
    \itemize{
      \item
      If \code{continuous=TRUE} (the default), smoothing is performed
      using the \dQuote{equal-split continuous} rule described in
      Section 9.2.3 of Okabe and Sugihara (2012).
      The resulting function is continuous on the linear network.
      \item
      If \code{continuous=FALSE}, smoothing is performed
      using the \dQuote{equal-split discontinuous} rule described in
      Section 9.2.2 of Okabe and Sugihara (2012). The
      resulting function is continuous except at the network vertices.
      \item
      In the default case
      (where \code{distance="path"} and
      \code{continuous=TRUE} and \code{kernel="gaussian"},
      computation is performed rapidly by solving the classical heat equation
      on the network, as described in McSwiggan et al (2016).
      The arguments are passed to \code{\link{densityHeat}} which performs
      the computation.
      Computational time is short, but increases quadratically
      with \code{sigma}.
      \item
      In all other cases, computation is performed by path-tracing
      as described in Okabe and Sugihara (2012);
      the arguments are passed to \code{\link{densityEqualSplit}}
      which performs the computation.
      Computation time can be extremely long, and
      increases exponentially with \code{sigma}.
    }
    \item If \code{distance="euclidean"}, the smoothing is performed
    using a two-dimensional kernel. The arguments are passed to
    \code{\link{densityQuick.lpp}} to perform the computation.
    Computation time is very short.
    See the help for \code{\link{densityQuick.lpp}} for further details.
  }
  There is also a method for split point patterns on a linear network
  (class \code{"splitppx"}) which will return a list of pixel images.
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
  \code{\link{lpp}},
  \code{\link{linim}},
  \code{\link{densityQuick.lpp}}
}
\examples{
  X <- runiflpp(3, simplenet)
  D <- density(X, 0.2, verbose=FALSE)
  plot(D, style="w", main="", adjust=2)
  Dq <- density(X, 0.2, distance="euclidean")
  plot(Dq, style="w", main="", adjust=2)
  Dw <- density(X, 0.2, weights=c(1,2,-1), verbose=FALSE)
  De <- density(X, 0.2, kernel="epanechnikov", verbose=FALSE)
  Ded <- density(X, 0.2, kernel="epanechnikov", continuous=FALSE, verbose=FALSE)
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
