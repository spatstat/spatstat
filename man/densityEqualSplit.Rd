\name{densityEqualSplit}
\alias{densityEqualSplit}
\title{
  Equal-Split Algorithm for Kernel Density on a Network
}
\description{
  Computes a kernel density estimate on a linear network
  using the Okabe-Sugihara equal-split algorithms.
}
\usage{
 densityEqualSplit(x, sigma = NULL, \dots, weights = NULL,
                   kernel = "epanechnikov", continuous = TRUE,
                   epsilon = 1e-06, verbose = TRUE, debug = FALSE, savehistory = TRUE)
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
    resolution of the result.
  }
  \item{weights}{
    Optional. Numeric vector of weights associated with the
    points of \code{x}. Weights may be positive, negative or zero.
  }
  \item{kernel}{
    Character string specifying the smoothing kernel.
    See \code{\link{dkernel}} for possible options.
  }
  \item{continuous}{
    Logical value indicating whether to compute the
    \dQuote{equal-split continuous} smoother (\code{continuous=TRUE}, the
    default) or the \dQuote{equal-split discontinuous} smoother
    (\code{continuous=FALSE}). 
  }
  \item{epsilon}{
    Tolerance value. A tail of the kernel with total mass
    less than \code{epsilon} may be deleted.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
  \item{debug}{
    Logical value indicating whether to print debugging information.
  }
  \item{savehistory}{
    Logical value indicating whether to save the entire history of the
    algorithm, for the purposes of evaluating performance.
  }
}
\details{
  Kernel smoothing is applied to the points of \code{x}
  using a kernel based on path distances in the network.
  The result is a pixel image on the linear network (class
  \code{"linim"}) which can be plotted.

  Smoothing is performed using one of the \dQuote{equal-split} rules described in
  Okabe and Sugihara (2012).
  
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
    resulting function is not continuous.
  }
  Computation is performed by path-tracing
  as described in Okabe and Sugihara (2012).

  It is advisable to choose a kernel with bounded support
  such as \code{kernel="epanechnikov"}.
  With a Gaussian kernel, computation time can be long, and
  increases exponentially with \code{sigma}.

  Faster algorithms are available through \code{\link{density.lpp}}.
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
  De <- density(X, 0.2, kernel="epanechnikov", verbose=FALSE)
  Ded <- density(X, 0.2, kernel="epanechnikov", continuous=FALSE, verbose=FALSE)
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
