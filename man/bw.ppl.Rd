\name{bw.ppl}
\alias{bw.ppl}
\title{
  Likelihood Cross Validation Bandwidth Selection for Kernel Density
}
\description{
  Uses likelihood cross-validation to select a smoothing bandwidth
  for the kernel estimation of point process intensity.
}
\usage{
   bw.ppl(X, ..., srange=NULL, ns=16)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{Ignored.}
  \item{srange}{
    Optional numeric vector of length 2 giving the
    range of values of bandwidth to be searched.
  }
  \item{ns}{
    Optional integer giving the number of values of
    bandwidth to search.
  }
}
\details{
  This function selects an appropriate bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.ppp}}.

  The bandwidth \eqn{\sigma}{sigma} is chosen to 
  maximise the point process likelihood cross-validation criterion
  \deqn{
    \mbox{LCV}(\sigma) =
    \sum_i \log\hat\lambda_{-i}(x_i) - \int_W \hat\lambda(u) \, {\rm d}u
  }{
    LCV(sigma) = sum[i] log(lambda[-i](x[i])) - integral[W] lambda(u) du
  }
  where the sum is taken over all the data points \eqn{x_i}{x[i]},
  where \eqn{\hat\lambda_{-i}(x_i)}{lambda[-i](x_i)} is the
  leave-one-out kernel-smoothing estimate of the intensity at
  \eqn{x_i}{x[i]} with smoothing bandwidth \eqn{\sigma}{sigma},
  and \eqn{\hat\lambda(u)}{lambda(u)} is the kernel-smoothing estimate
  of the intensity at a spatial location \eqn{u} with smoothing
  bandwidth \eqn{\sigma}{sigma}.
  See Loader(1999, Section 5.3).

  The value of \eqn{\mbox{LCV}(\sigma)}{LCV(sigma)} is computed
  directly, using \code{\link{density.ppp}}, 
  for \code{ns} different values of \eqn{\sigma}{sigma}
  between \code{srange[1]} and \code{srange[2]}.

  The result is a numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted to show the (rescaled) mean-square error
  as a function of \code{sigma}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{bw.diggle}},
  \code{\link{bw.scott}}
}
\examples{
  \donttest{
    b <- bw.ppl(redwood)
    plot(b, main="Likelihood cross validation for redwoods")
    plot(density(redwood, b))
  }
  \testonly{
    b <- bw.ppl(redwood, srange=c(0.03, 0.07), ns=2)
  }
}
\references{
  Loader, C. (1999)
  \emph{Local Regression and Likelihood}.
  Springer, New York.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
