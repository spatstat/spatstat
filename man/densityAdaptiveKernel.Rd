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
    Point pattern (object of class \code{"ppp"}).
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
    \code{\link{bw.abram}} to compute the smoothing bandwidths
    if \code{bw} is missing,
    or passed to \code{\link[spatstat]{as.mask}}
    to control the spatial resolution of the result.
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
  This function computes a spatially-adaptive kernel estimate of the
  spatially-varying intensity from the point pattern \code{X}
  using the partitioning technique of Davies and Baddeley (2018).

  The argument \code{bw} specifies the smoothing bandwidths to be
  applied to each of the points in \code{X}. It may be a numeric vector
  of bandwidth values, or a pixel image or function yielding the
  bandwidth values.

  If the points of \code{X} are \eqn{x_1,\ldots,x_n}{x[1], ..., x[n]}
  and the corresponding bandwidths are
  \eqn{\sigma_1,\ldots,\sigma_n}{\sigma[1], ..., \sigma[n]}
  then the adaptive kernel estimate of intensity at a location \eqn{u} is
  \deqn{
    \hat\lambda(u) = \sum_{i=1}^n k(u, x_i, \sigma_i)
  }{
    \lambda(u) = sum[i] e(x[i], k(u, x[i], \sigma[i])
  }
  where \eqn{k(u, v, \sigma)} is the value at \eqn{u}
  of the (possibly edge-corrected) smoothing kernel with bandwidth \eqn{\sigma}
  induced by a data point at \eqn{v}.

  Exact computation of the estimate above can be time-consuming:
  it takes \eqn{n} times longer than fixed-bandwidth smoothing.
  
  The partitioning method of Davies and Baddeley (2018)
  accelerates this computation by partitioning the range of
  bandwidths into \code{ngroups} intervals,
  correspondingly subdividing the points of the pattern \code{X} into
  \code{ngroups} sub-patterns according to bandwidth,
  and applying fixed-bandwidth smoothing to each sub-pattern.

  The default value of \code{ngroups} is the integer part of the square root of
  the number of points in \code{X}, so that the computation time is
  only about \eqn{\sqrt{n}}{sqrt(n)} times slower than fixed-bandwidth
  smoothing. Any positive value of \code{ngroups}
  can be specified by the user. Specifying \code{ngroups=Inf} enforces exact
  computation of the estimate without partitioning. Specifying
  \code{ngroups=1} is the same as fixed-bandwidth smoothing with
  bandwidth \code{sigma=median(bw)}.
}
\section{Bandwidths and Bandwidth Selection}{
  The function \code{densityAdaptiveKernel}
  computes one adaptive estimate of the intensity,
  determined by the smoothing bandwidth values \code{bw}.

  Typically the bandwidth values are computed by first computing
  a pilot estimate of the intensity, then using \code{\link{bw.abram}}
  to compute the vector of bandwidths according to Abramson's rule.
  This involves specifying a global bandwidth \code{h0}.
  
  The default bandwidths may work well in many contexts, but for optimal
  bandwidth selection, this calculation should be performed repeatedly with
  different values of \code{h0} to optimise the value of \code{h0}.
  This can be computationally demanding; we recommend
  the function \code{multiscale.density} in the \pkg{sparr} package
  which supports much faster bandwidth selection, using the FFT
  method of Davies and Baddeley (2018).
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

  See the function \code{bivariate.density} in the \pkg{sparr} package
  for a more flexible implementation, and
  \code{multiscale.density} for an implementation that is more
  efficient for bandwidth selection.
}
\examples{
  Z <- densityAdaptiveKernel(redwood, h0=0.1)
  plot(Z, main="Adaptive kernel estimate")
  points(redwood, col="white")
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
