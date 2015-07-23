\name{dkernel}
\alias{dkernel}
\alias{pkernel}
\alias{qkernel}
\alias{rkernel}
\title{Kernel distributions and random generation}
\description{Density, distribution function, quantile function and random
  generation for several distributions used in kernel estimation
  for numerical data.
}
\usage{
dkernel(x, kernel = "gaussian", mean = 0, sd = 1)
pkernel(q, kernel = "gaussian", mean = 0, sd = 1, lower.tail = TRUE)
qkernel(p, kernel = "gaussian", mean = 0, sd = 1, lower.tail = TRUE)
rkernel(n, kernel = "gaussian", mean = 0, sd = 1)
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{kernel}{
    String name of the kernel.
    Options are
    \code{"gaussian"}, \code{"rectangular"},
    \code{"triangular"},
    \code{"epanechnikov"},
    \code{"biweight"},
    \code{"cosine"} and \code{"optcosine"}.
    (Partial matching is used).
  }
  \item{n}{Number of observations.}
  \item{mean}{Mean of distribution.}
  \item{sd}{Standard deviation of distribution.}
  \item{lower.tail}{logical; if \code{TRUE} (the default),
    then probabilities are \eqn{P(X \le x)}{P[X <= x]},
    otherwise, \eqn{P(X > x)}.
  }
}
\details{
  These functions give the
  probability density, cumulative distribution function,
  quantile function and random generation for several
  distributions used in kernel estimation for one-dimensional
  (numerical) data.

  The available kernels are those used in \code{\link[stats]{density.default}},
  namely \code{"gaussian"}, \code{"rectangular"},
    \code{"triangular"},
    \code{"epanechnikov"},
    \code{"biweight"},
    \code{"cosine"} and \code{"optcosine"}.
    For more information about these kernels,
    see \code{\link[stats]{density.default}}. 
  
  \code{dkernel} gives the probability density,
  \code{pkernel} gives the cumulative distribution function,
  \code{qkernel} gives the quantile function,
  and \code{rkernel} generates random deviates.
}
\value{
  A numeric vector.
  For \code{dkernel}, a vector of the same length as \code{x}
  containing the corresponding values of the probability density.
  For \code{pkernel}, a vector of the same length as \code{x}
  containing the corresponding values of the cumulative distribution function.
  For \code{qkernel}, a vector of the same length as \code{p}
  containing the corresponding quantiles.
  For \code{rkernel}, a vector of length \code{n}
  containing randomly generated values.
}
\examples{
  x <- seq(-3,3,length=100)
  plot(x, dkernel(x, "epa"), type="l",
           main=c("Epanechnikov kernel", "probability density"))
  plot(x, pkernel(x, "opt"), type="l",
           main=c("OptCosine kernel", "cumulative distribution function"))
  p <- seq(0,1, length=256)
  plot(p, qkernel(p, "biw"), type="l",
           main=c("Biweight kernel", "cumulative distribution function"))
  y <- rkernel(100, "tri")
  hist(y, main="Random variates from triangular density")
  rug(y)
}
\seealso{
  \code{\link[stats]{density.default}},
  \code{\link{kernel.factor}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Martin Hazelton
}
\keyword{methods}
\keyword{nonparametric}
\keyword{smooth}
