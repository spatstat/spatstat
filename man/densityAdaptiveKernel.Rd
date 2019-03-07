\name{densityAdaptiveKernel}
\alias{densityAdaptiveKernel}
\alias{densityAdaptiveKernel.ppp}
\title{Adaptive Kernel Estimate of Intensity of Point Pattern}
\description{
  Computes an adaptive estimate of the intensity function of a point
  pattern using a variable-bandwidth smoothing kernel.
}
\usage{
densityAdaptiveKernel(X, \dots)

\method{densityAdaptiveKernel}{ppp}(X, bw, \dots,
       weights=NULL,       
       at=c("pixels", "points"),
       edge=TRUE, ngroups) 
}
\arguments{
  \item{X}{
    Point pattern in two dimensions.
  }
  \item{bw}{
    Numeric vector of smoothing bandwidths for each point in \code{X},
    or a pixel image giving the smoothing bandwidth at each spatial
    location, or a spatial function of class \code{"funxy"} giving the
    smoothing bandwidth at each location.
    The default is to compute bandwidths using \code{\link{bw.abram}}.
  }
  \item{\dots}{
    Arguments passed to
    \code{\link{bw.abram}} to compute the smoothing bandwidths,
    or passed to \code{\link[spatstat]{as.mask}}
    to control the spatial resolution.
  }
  \item{weights}{
    Optional vector of numeric weights for the points of \code{X}.
  }
  \item{at}{
    String specifying whether to compute the intensity values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{x} (\code{at="points"}).
  }
  \item{edge}{
    Logical value indicating whether to perform edge correction.
  }
  \item{ngroups}{
    Number of groups into which the bandwidth values should be
    partitioned and discretised.
  }
}
\details{
  To be written.

  Uses the splitting technique of Davies and Baddeley (2018).
}
\value{
  If \code{at="pixels"} (the default), the result is a pixel image.
  If \code{at="points"}, the result is a numeric vector with one entry
  for each data point in \code{X}.
}
\author{
  \adrian and Tilman Davies.
}
\references{
  Davies, T.M. and Baddeley, A. (2018)
  Fast computation of spatially adaptive kernel estimates.
  \emph{Statistics and Computing}, \bold{28}(4), 937-956.\cr

  Hall, P. and Marron, J.S. (1988)
  Variable window width kernel density estimates of probability
  densities.
  \emph{Probability Theory and Related Fields}, \bold{80}, 37-49.\cr

  Silverman, B.W. (1986)
  \emph{Density Estimation for Statistics and Data Analysis}.
  Chapman and Hall, New York.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{adaptive.density}},
  \code{\link{densityVoronoi}},
  \code{\link{im.object}}.
}
\examples{
  Z <- densityAdaptiveKernel(redwood, h0=0.1)
  plot(Z, main="Adaptive kernel estimate")
  points(redwood, col="white")
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
