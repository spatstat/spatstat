\name{density.psp}
\alias{density.psp}
\alias{density.linnet}
\title{Kernel Smoothing of Line Segment Pattern or Linear Network}
\description{
  Compute a kernel smoothed intensity function from a line segment
  pattern or a linear network.
}
\usage{
  \method{density}{psp}(x, sigma, \dots, edge=TRUE,
                   method=c("FFT", "C", "interpreted"),
                   at=NULL)

  \method{density}{linnet}(x, \dots)
}
\arguments{
  \item{x}{
    Line segment pattern (object of class \code{"psp"})
    or linear network (object of class \code{"linnet"})
    to be smoothed.
  }
  \item{sigma}{
    Standard deviation of isotropic Gaussian smoothing kernel.
  }
  \item{\dots}{
    Extra arguments, including arguments passed to \code{\link{as.mask}}
    to determine the resolution of the resulting image.
  }
  \item{edge}{
    Logical flag indicating whether to apply edge correction.
  }
  \item{method}{
    Character string (partially matched) specifying the method of
    computation. Option \code{"FFT"} is the fastest, while
    \code{"C"} is the most accurate.
  }
  \item{at}{
    Optional. An object specifying the locations where density values
    should be computed. Either a window (object of class \code{"owin"})
    or a point pattern (object of class \code{"ppp"} or \code{"lpp"}).
  }
}
\value{
  A pixel image (object of class \code{"im"})
  or a numeric vector.
}
\details{
  These are methods for the generic function \code{\link{density}}
  for the classes \code{"psp"} (line segment patterns)
  and \code{"linnet"} (linear networks).
  If \code{x} is a linear network, it is first converted
  to a line segment pattern.

  A kernel estimate of the intensity of the line segment pattern
  is computed. The result is 
  the convolution of the isotropic Gaussian kernel, of
  standard deviation \code{sigma}, with the line segments.
  The result is computed as follows:
  \itemize{
    \item if \code{method="FFT"} (the default),
    the line segments are discretised
    using \code{\link{pixellate.psp}}, then the Fast Fourier Transform
    is used to calculate the convolution. This method is the fastest,
    but is slightly less accurate. Accuracy can be improved by
    increasing pixel resolution.
    \item if \code{method="C"} the exact value of the convolution at the
    centre of each pixel is computed analytically using \code{C} code;
    \item if \code{method="interpreted"},
    the exact value of the convolution at the
    centre of each pixel is computed analytically using \code{R} code.
    This method is the slowest.
  }
  If \code{edge=TRUE} this result is adjusted for edge effects
  by dividing it by the convolution of the same Gaussian kernel
  with the observation window.

  If the argument \code{at} is given, then it specifies the locations
  where density values should be computed.
  \itemize{
    \item If \code{at} is a window,
    then the window is converted to a binary mask
    using the arguments \code{\dots}, and density values are computed at
    the centre of each pixel in this mask. The result is a pixel image.
    \item
    If \code{at} is a point pattern, then density values are computed at
    each point location, and the result is a numeric vector.
  }
}
\seealso{
  \code{\link{psp.object}},
  \code{\link{im.object}},
  \code{\link{density}}
}
\examples{
  L <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
  D <- density(L, sigma=0.03)
  plot(D, main="density(L)")
  plot(L, add=TRUE)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
