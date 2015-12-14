\name{Kest.fft}
\alias{Kest.fft}
\title{K-function using FFT}
\description{
Estimates the reduced second moment function \eqn{K(r)} 
from a point pattern in a window of arbitrary shape,
using the Fast Fourier Transform.
}
\usage{
  Kest.fft(X, sigma, r=NULL, \dots, breaks=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{K(r)} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{sigma}{
    Standard deviation of the isotropic Gaussian
    smoothing kernel.
  }
  \item{r}{
    Optional. Vector of values for the argument \eqn{r} at which \eqn{K(r)} 
    should be evaluated. There is a sensible default.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} determining the
    spatial resolution for the FFT calculation.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
  
  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{K} has been  estimated
  }
  \item{border}{the estimates of \eqn{K(r)} for these values of \eqn{r}
  }
  \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
    for a stationary Poisson process
  }
}
\details{
  This is an alternative to the function \code{\link{Kest}}
  for estimating the \eqn{K} function. It may be useful for
  very large patterns of points.

  Whereas \code{\link{Kest}} computes the distance between
  each pair of points analytically, this function discretises the
  point pattern onto a rectangular pixel raster and applies
  Fast Fourier Transform techniques to estimate \eqn{K(t)}.
  The hard work is done by the function \code{\link{Kmeasure}}.

  The result is an approximation whose accuracy depends on the
  resolution of the pixel raster. The resolution is controlled
  by the arguments \code{\dots}, or by setting the parameter \code{npixel} in 
  \code{\link{spatstat.options}}.
}
\references{
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Academic Press, 1983.

  Ohser, J. (1983)
  On estimators for the reduced second moment measure of
  point processes. \emph{Mathematische Operationsforschung und
  Statistik, series Statistics}, \bold{14}, 63 -- 71.
    
  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J. (1995)
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
} 
\seealso{
  \code{\link{Kest}},
  \code{\link{Kmeasure}},
  \code{\link{spatstat.options}}
}
\examples{
 pp <- runifpoint(10000)
 \dontrun{
 spatstat.options(npixel=512)
 }
 \testonly{
  op <- spatstat.options(npixel=125)
 }
 Kpp <- Kest.fft(pp, 0.01)
 plot(Kpp)
 \testonly{spatstat.options(op)}
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{nonparametric}
 
 
