\name{density.lpp}
\alias{density.lpp}
\alias{density.splitppx}
\title{
  Kernel Estimate of Intensity on a Linear Network
}
\description{
  Estimates the intensity of a point process on a linear network
  by applying kernel smoothing to the point pattern data,
  using the equal-split continuous algorithm.
}
\usage{
\method{density}{lpp}(x, sigma, \dots,
        kernel="gaussian", continuous=TRUE,
        epsilon = 1e-06, verbose = TRUE, debug = FALSE, savehistory = TRUE)

\method{density}{splitppx}(x, sigma, \dots)
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
  using one of the algorithms described in Okabe and Sugihara (2012).
  The result is a pixel image on the linear network (class
  \code{"linim"}) which can be plotted.

  If \code{continuous=TRUE} (the default), smoothing is performed
  using the \dQuote{equal-split continuous} rule described in
  Section 9.2.3 of Okabe and Sugihara (2012). The resulting function
  is continuous on the linear network.

  If \code{continuous=FALSE}, smoothing is performed
  using the \dQuote{equal-split discontinuous} rule described in
  Section 9.2.2 of Okabe and Sugihara (2012). This algorithm is
  faster, but the resulting function is not continuous.  

  There is also a method for split point patterns on a linear network
  (class \code{"splitppx"}) which will return a list of pixel images.
}
\section{WARNING}{
  \bold{THIS ALGORITHM CAN BE EXTREMELY SLOW} for large values of \code{sigma}.

  The computational complexity increases exponentially with
  \code{sigma}.
  
  You Have Been Warned.
}
\value{
  Pixel image on the linear network (class
  \code{"linim"}).
}
\references{
  Okabe, A. and Sugihara, K. (2012)
  \emph{Spatial analysis along networks}.
  Wiley.
}
\author{
  Adrian Baddeley and Greg McSwiggan.
}
\seealso{
  \code{\link{lpp}},
  \code{\link{linim}}
}
\examples{
  X <- runiflpp(3, simplenet)
  D <- density(X, 0.2, verbose=FALSE)
  plot(D, style="w", main="", adjust=2)
  De <- density(X, 0.2, kernel="epanechnikov", verbose=FALSE)
  Ded <- density(X, 0.2, kernel="epanechnikov", continuous=FALSE, verbose=FALSE)
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
