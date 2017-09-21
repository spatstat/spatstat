\name{density.psp}
\alias{density.psp}
\title{Kernel Smoothing of Line Segment Pattern}
\description{
  Compute a kernel smoothed intensity function from a line segment pattern.
}
\usage{
  \method{density}{psp}(x, sigma, \dots, edge=TRUE,
                   method=c("FFT", "C", "interpreted"))
}
\arguments{
  \item{x}{
    Line segment pattern (object of class \code{"psp"}) to be smoothed.
  }
  \item{sigma}{
    Standard deviation of isotropic Gaussian smoothing kernel.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{as.mask}} which determine
    the resolution of the resulting image.
  }
  \item{edge}{
    Logical flag indicating whether to apply edge correction.
  }
  \item{method}{
    Character string (partially matched) specifying the method of
    computation. Option \code{"FFT"} is the fastest, while
    \code{"C"} is the most accurate.
  }
}
\value{
  A pixel image (object of class \code{"im"}).
}
\details{
  This is a method for the generic function \code{\link{density}}.

  A kernel estimate of the intensity of the line segment pattern
  is computed. The result is 
  the convolution of the isotropic Gaussian kernel, of
  standard deviation \code{sigma}, with the line segments.
  The result is computed as follows:
  \itemize{
    \item if \code{method="FFT"}, the line segments are discretised
    using \code{\link{pixellate.psp}}, then the Fast Fourier Transform
    is used to calculate the convolution. This method is the fastest,
    but is slightly less accurate.
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
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
