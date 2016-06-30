\name{bw.scott}
\alias{bw.scott}
\title{
  Scott's Rule for Bandwidth Selection for Kernel Density
}
\description{
  Use Scott's rule of thumb to determine the smoothing bandwidth
  for the kernel estimation of point process intensity.
}
\usage{
   bw.scott(X)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
}
\details{
  This function selects a bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.ppp}}.

  The bandwidth \eqn{\sigma}{sigma} is computed by the rule of thumb
  of Scott (1992, page 152). It is very fast to compute.

  This rule is designed for density
  estimation, and typically produces a larger bandwidth
  than \code{\link{bw.diggle}}. It is useful for estimating
  gradual trend. 
}
\value{
  A numerical vector of two elements giving the selected
  bandwidths in the \code{x} and \code{y} directions.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{bw.diggle}},
  \code{\link{bw.ppl}},
  \code{\link{bw.frac}}.
}
\examples{
  data(lansing)
  attach(split(lansing))
  b <- bw.scott(hickory)
  b
  \donttest{
   plot(density(hickory, b))
  }
}
\references{
  Scott, D.W. (1992)
  \emph{Multivariate Density Estimation. Theory, Practice and
    Visualization}. 
  New York: Wiley.
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
