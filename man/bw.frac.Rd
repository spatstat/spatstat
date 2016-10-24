\name{bw.frac}
\alias{bw.frac}
\title{
  Bandwidth Selection Based on Window Geometry
}
\description{
  Select a smoothing bandwidth for smoothing a point pattern,
  based only on the geometry of the spatial window.
  The bandwidth is a specified quantile of the distance
  between two independent random points in the window.
}
\usage{
   bw.frac(X, \dots, f=1/4)
}
\arguments{
  \item{X}{
    A window (object of class \code{"owin"}) or
    point pattern (object of class \code{"ppp"})
    or other data which can be converted to a window
    using \code{\link{as.owin}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{distcdf}}.
  }
  \item{f}{
    Probability value (between 0 and 1)
    determining the quantile of the distribution.
  }
}
\details{
  This function selects an appropriate bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.ppp}}.

  The bandwidth \eqn{\sigma}{\sigma} is computed as a 
  quantile of the distance between two independent random points
  in the window. The default is the lower quartile of this
  distribution.

  If \eqn{F(r)} is the cumulative distribution function of the
  distance between two independent random points uniformly distributed
  in the window, then the value returned is the quantile
  with probability \eqn{f}. That is, the bandwidth is 
  the value \eqn{r} such that \eqn{F(r) = f}.

  The cumulative distribution function  \eqn{F(r)} is
  computed using \code{\link{distcdf}}. We then
  we compute the smallest number \eqn{r}
  such that \eqn{F(r) \ge f}{F(r) >= f}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.frac"}
  which can be plotted to show the cumulative distribution function
  and the selected quantile.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{bw.diggle}},
  \code{\link{bw.ppl}},
  \code{\link{bw.relrisk}},
  \code{\link{bw.scott}},
  \code{\link{bw.smoothppp}},
  \code{\link{bw.stoyan}}
}
\examples{
  h <- bw.frac(letterR)
  h
  plot(h, main="bw.frac(letterR)")
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
