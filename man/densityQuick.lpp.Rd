\name{densityQuick.lpp}
\alias{densityQuick.lpp}
\title{
  Kernel Estimation of Intensity on a Network using a 2D Kernel
}
\description{
  Estimates the intensity of a point process on a linear network
  using a two-dimensional smoothing kernel.
}
\usage{
densityQuick.lpp(X, sigma=NULL, \dots,
         kernel="gaussian",
         at = c("pixels", "points"),
         what = c("estimate", "se", "var"),
         leaveoneout = TRUE,
         diggle = FALSE,
         edge2D = FALSE,
         weights = NULL,
         positive = FALSE)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{sigma}{
    Smoothing bandwidth. A single numeric value,
    in the same units as the coordinates of \code{X}.
    Alternatively \code{sigma} may be a function which selects a
    bandwidth when applied to \code{X},
    for example, \code{\link{bw.scott.iso}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{as.mask}}
    to determine the pixel resolution.
  }
  \item{kernel}{
    String (partially matched) specifying the smoothing kernel.
    Current options are \code{"gaussian"}, \code{"epanechnikov"},
    \code{"quartic"} or \code{"disc"}.
  }
  \item{at}{
    String (partially matched)
    specifying whether to compute the intensity values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{x} (\code{at="points"}).
  }
  \item{what}{
    String (partially matched)
    specifying whether to calculate the intensity estimate,
    or its estimated standard error, or its estimated variance.
  }
  \item{leaveoneout}{
    Logical value indicating whether to compute a leave-one-out
    estimator. Applicable only when \code{at="points"}.
  }
  \item{diggle}{
    Logical value specifying whether to use the \sQuote{Diggle}
    correction. 
  }
  \item{edge2D}{
    Logical value specifying whether to apply the usual
    two-dimensional edge correction procedure to the numerator and
    denominator of the estimate.
}
  \item{weights}{
    Optional weights to be attached to the points.
    A numeric vector, an \code{expression}, or a pixel image.
  }
  \item{positive}{
    Logical value indicating whether to force the resulting values to be
    positive. Default is \code{FALSE} for the sake of speed.
  }
}
\details{
  Kernel smoothing is applied to the points of \code{x}
  using a two-dimensional Gaussian kernel, as described in Rakshit et al (2019).
  The result is a pixel image on the linear network (class
  \code{"linim"}) which can be plotted.

  Other techniques for kernel smoothing on a network are implemented in
  \code{\link{density.lpp}}. The main advantages of using a
  two-dimensional kernel are very fast computation and
  insensitivity to changes in the network geometry. The main
  disadvantage is that it ignores the connectivity of the network.
  See Rakshit et al (2019) for further explanation.
}
\value{
  Image on a linear network (object of class \code{"linim"}).
}
\references{
  Rakshit, S., Davies, T., Moradi, M., 
  McSwiggan, G., Nair, G., Mateu, J. and Baddeley, A. (2019)
  Fast kernel smoothing of point patterns on a large
  network using 2D convolution.
  \emph{International Statistical Review}.
  In press. Published online 06 June 2019.
  DOI: 10.1111/insr.12327.
}
\author{
  Adrian Baddeley, Suman Rakshit and Tilman Davies
}
\seealso{
  \code{\link{density.lpp}}, the main function for density estimation
  on a network.

  \code{\link{bw.scott}}, \code{\link{bw.scott.iso}} for bandwidth selection.
}
\examples{
   X <- unmark(chicago)
   plot(densityQuick.lpp(X, 500))
   plot(densityQuick.lpp(X, 500, diggle=TRUE))
   plot(densityQuick.lpp(X, bw.scott.iso))
   plot(densityQuick.lpp(X, 500, what="se"))
}
\keyword{spatial}
\keyword{nonparametric}

