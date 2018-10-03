\name{bw.CvL}
\alias{bw.CvL}
\title{
  Cronie and van Lieshout's Criterion for Bandwidth Selection for Kernel Density
}
\description{
  Uses Cronie and van Lieshout's criterion based on Cambell's formula to select
  a smoothing bandwidth for the kernel estimation of point process intensity.
}
\usage{
   bw.CvL(X, \dots, srange = NULL, ns = 16, sigma = NULL)
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
  \item{sigma}{
    Optional. Vector of values of the bandwidth to be searched.
    Overrides the values of \code{ns} and \code{srange}.
  }
}
\details{
  This function selects an appropriate bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.ppp}}.

  The bandwidth \eqn{\sigma}{\sigma} is chosen to 
  minimise the discrepancy between the area of the observation window and the
  sum of reciprocal estimated intensity values at the points of the point process
  \deqn{
    \mbox{CvL}(\sigma) =
    (|W| - \sum_i 1/\hat\lambda(x_i))^2
  }{
    CvL(\sigma) = (|W| - sum[i] 1/\lambda(x[i]))^2
  }
  where the sum is taken over all the data points \eqn{x_i}{x[i]},
  and where \eqn{\hat\lambda(x_i)}{\lambda(x[i])} is the
  kernel-smoothing estimate of the intensity at
  \eqn{x_i}{x[i]} with smoothing bandwidth \eqn{\sigma}{\sigma}.

  The value of \eqn{\mbox{CvL}(\sigma)}{CvL(\sigma)} is computed
  directly, using \code{\link{density.ppp}}, 
  for \code{ns} different values of \eqn{\sigma}{\sigma}
  between \code{srange[1]} and \code{srange[2]}.

  The result is a numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted to show the bandwidth selection criterion
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
  \code{\link{bw.scott}},
  \code{\link{bw.ppl}},
  \code{\link{bw.frac}}.
}
\examples{
  if(interactive()) {
    b <- bw.CvL(redwood)
    b
    plot(b, main="Cronie and van Lieshout bandwidth criterion for redwoods")
    plot(density(redwood, b))
    plot(density(redwood, bw.CvL))
  }
  \testonly{
    b <- bw.CvL(redwood, srange=c(0.03, 0.07), ns=2)
  }
}
\references{
 Cronie, O and Van Lieshout, M N M (2018) A non-model-based approach to
 bandwidth selection for kernel estimators of spatial intensity functions,
 \emph{Biometrika}, \bold{105}, 455-462.
}
\author{Ottmar Cronie
  \email{ottmar.cronie@umu.se}
  and Marie-Colette van Lieshout
  \email{Marie-Colette.van.Lieshout@cwi.nl}
  adapted for \pkg{spatstat} by
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
