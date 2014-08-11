\name{bw.diggle}
\alias{bw.diggle}
\title{
  Cross Validated Bandwidth Selection for Kernel Density
}
\description{
  Uses cross-validation to select a smoothing bandwidth
  for the kernel estimation of point process intensity.
}
\usage{
   bw.diggle(X, ...)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{Ignored.}
}
\details{
  This function selects an appropriate bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.ppp}}.

  The bandwidth \eqn{\sigma}{sigma} is chosen to 
  minimise the mean-square error criterion defined by Diggle (1985).
  The algorithm computes the mean-square error by 
  the method of Berman and Diggle (1989).
  See Diggle (2003, pages 115-118) for a summary of this method.

  The result is a numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted to show the (rescaled) mean-square error
  as a function of \code{sigma}.
}
\section{Definition of bandwidth}{
  The smoothing parameter \code{sigma} returned by \code{bw.diggle}
  (and displayed on the horizontal axis of the plot)
  corresponds to \code{h/2}, where \code{h} is the smoothing
  parameter described in Diggle (2003, pages 116-118) and
  Berman and Diggle (1989).
  In those references, the smoothing kernel 
  is the uniform density on the disc of radius \code{h}. In
  \code{\link{density.ppp}}, the smoothing kernel is the
  isotropic Gaussian density with standard deviation \code{sigma}.
  When replacing one kernel by another, the usual
  practice is to adjust the bandwidths so that the kernels have equal
  variance (cf. Diggle 2003, page 118). This implies that \code{sigma = h/2}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{bw.ppl}},
  \code{\link{bw.scott}}
}
\examples{
  data(lansing)
  attach(split(lansing))
  b <- bw.diggle(hickory)
  plot(b, ylim=c(-2, 0), main="Cross validation for hickories")
  \donttest{
   plot(density(hickory, b))
  }
}
\references{
  Berman, M. and Diggle, P. (1989)
  Estimating weighted integrals of the
  second-order intensity of a spatial point process.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{51}, 81--92.

  Diggle, P.J. (1985)
  A kernel method for smoothing point process data.
  \emph{Applied Statistics} (Journal of the Royal Statistical Society,
  Series C) \bold{34} (1985) 138--147.

  Diggle, P.J. (2003)
  \emph{Statistical analysis of spatial point patterns},
  Second edition. Arnold.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
