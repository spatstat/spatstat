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
       bwrange = range(bw), 
       edge=TRUE, diggle=TRUE, weights=NULL, ngroups=NULL,
       nbw = 64, what = c("smooth", "all"),
       ncoarse=NULL, taperfac=c(5,10), verbose=TRUE)
}
\arguments{
  \item{X}{
    Point pattern in two dimensions.
  }
  \item{bw}{
    Numeric vector of smoothing bandwidths for each point in \code{X},
    or a pixel image giving the smoothing bandwidth at each spatial
    location.
    If missing, a bandwidth image is computed using
    \code{\link{bw.abram}}.
  }
  \item{bwrange}{
    Range of possible bandwidths. A numeric vector of length 2
    giving the minimum and maximum possible bandwidths.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{as.mask}}
    controlling the spatial resolution, or passed to
    \code{\link{bw.abram}} to compute the smoothing bandwidths.
  }
  \item{weights}{
    Optional vector of numeric weights for the points of \code{X}.
  }
  \item{ngroups}{
    Optional.
    Number of groups into which the bandwidth values should be partitioned.
    If \code{ngroups = NULL} (the default) the bandwidths are not partitioned,
    and the entire calculation is performed in one step.
    If \code{ngroups} is given, the bandwidths are partitioned into
    quantiles and the calculation is performed separately for each
    quantile group.
  }
  \item{nbw}{
    Integer. Total number of different bandwidth values 
    in the discretised grid.
  }
  \item{what}{
    Character string specifying whether to return only the final result,
    the variable-bandwidth kernel estimate of intensity
    (\code{what="smooth"}, the default)
    or a list containing many interim results of the calculation
    as well as this final result (\code{what="all"}).
  }
  \item{edge}{
    Logical value indicating whether to perform edge correction.
  }
  \item{diggle}{
    Logical value indicating whether to use the analogue of
    Diggle's edge correction (\code{diggle=TRUE}, the default) or the
    uniform edge correction (\code{diggle=FALSE}).
  }
  \item{ncoarse}{
    Numeric vector of length 3 providing the dimensions of
    a coarser 3D grid upon which the edge-correction factors
    should be computed, in the order of \code{c(nx,ny,nh)}.
  }
  \item{taperfac}{
    Positive numeric vector of length 2 providing control factors for
    the bandwidth tapering function; \code{[1]} for the lower taper,
    \code{[2]} for the upper taper.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  To be written.
}
\value{
  If \code{what="smooth"} (the default), the result is a pixel image.
  If \code{what="all"}, the result is a named list containing interim
  results of the calculation.
}
\author{
  Tilman Davies and \adrian.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{adaptive.density}},
  \code{\link{densityVoronoi}},
  \code{\link{im.object}}.
}
\examples{
  h0 <- bw.abram(redwood,h0=0.1,hp=0.07, at="pixels")
  nbw <- if(interactive()) 64 else 16
  Z <- densityAdaptiveKernel(redwood, h0, nbw=nbw)
  plot(Z, main="Adaptive kernel estimate")
  points(redwood, col="white")
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
