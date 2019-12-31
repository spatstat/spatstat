\name{bw.lppl}
\alias{bw.lppl}
\title{
  Likelihood Cross Validation Bandwidth Selection for Kernel Density on
  a Linear Network
}
\description{
  Uses likelihood cross-validation to select a smoothing bandwidth
  for the kernel estimation of point process intensity on a linear network.
}
\usage{
   bw.lppl(X, \dots, srange=NULL, ns=16, sigma=NULL, weights=NULL,
           distance="euclidean", shortcut=FALSE, warn=TRUE)
}
\arguments{
  \item{X}{
    A point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{srange}{
    Optional numeric vector of length 2 giving the
    range of values of bandwidth to be searched.
  }
  \item{ns}{
    Optional integer giving the number of values of
    bandwidth to search.
  }
  \item{sigma}{
    Optional. Vector of values of the bandwidth to be searched.
    Overrides the values of \code{ns} and \code{srange}.
  }
  \item{weights}{
    Optional. Numeric vector of weights for the points of \code{X}.
    Argument passed to \code{\link{density.lpp}}.
  }
  \item{distance}{
    Argument passed to \code{\link{density.lpp}} controlling the
    type of kernel estimator. 
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{density.lpp}}.
  }
  \item{shortcut}{
    Logical value indicating whether to speed up the calculation
    by omitting the integral term in the cross-validation criterion.
  }
  \item{warn}{
    Logical. If \code{TRUE}, issue a warning if the maximum of
    the cross-validation criterion occurs at one of the ends of the
    search interval.
  }
}
\details{
  This function selects an appropriate bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.lpp}}.

  The argument \code{X} should be a point pattern on a linear network
  (class \code{"lpp"}).

  The bandwidth \eqn{\sigma}{\sigma} is chosen to 
  maximise the point process likelihood cross-validation criterion
  \deqn{
    \mbox{LCV}(\sigma) =
    \sum_i \log\hat\lambda_{-i}(x_i) - \int_L \hat\lambda(u) \, {\rm d}u
  }{
    LCV(\sigma) = sum[i] log(\lambda[-i](x[i])) - integral[L] \lambda(u) du
  }
  where the sum is taken over all the data points \eqn{x_i}{x[i]},
  where \eqn{\hat\lambda_{-i}(x_i)}{\lambda[-i](x_i)} is the
  leave-one-out kernel-smoothing estimate of the intensity at
  \eqn{x_i}{x[i]} with smoothing bandwidth \eqn{\sigma}{\sigma},
  and \eqn{\hat\lambda(u)}{\lambda(u)} is the kernel-smoothing estimate
  of the intensity at a spatial location \eqn{u} with smoothing
  bandwidth \eqn{\sigma}{\sigma}.
  See Loader(1999, Section 5.3).

  The value of \eqn{\mbox{LCV}(\sigma)}{LCV(\sigma)} is computed
  directly, using \code{\link{density.lpp}}, 
  for \code{ns} different values of \eqn{\sigma}{\sigma}
  between \code{srange[1]} and \code{srange[2]}.

  The result is a numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted to show the (rescaled) mean-square error
  as a function of \code{sigma}.

  If \code{shortcut=TRUE}, the computation is accelerated by
  omitting the integral term in the equation above. This is valid
  because the integral is approximately constant.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\seealso{
  \code{\link{density.lpp}},
  \code{\link{bw.scott}}.

  For point patterns in two-dimensional space, use \code{\link{bw.ppl}}.
}
\examples{
  if(interactive()) {
    b <- bw.lppl(spiders)
    plot(b, main="Likelihood cross validation for spiders")
    plot(density(spiders, b, distance="e"))
  } else {
    b1 <- bw.lppl(spiders, ns=2)
    b2 <- bw.lppl(spiders, ns=2, shortcut=TRUE)
  }
}
\references{
  Loader, C. (1999)
  \emph{Local Regression and Likelihood}.
  Springer, New York.
}
\author{
  Suman Rakshit.
  Tweaks by \adrian.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
