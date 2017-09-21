\name{dmixpois}
\alias{dmixpois}
\alias{pmixpois}
\alias{qmixpois}
\alias{rmixpois}
\title{
  Mixed Poisson Distribution
}
\description{
  Density, distribution function, quantile function and random
  generation for a mixture of Poisson distributions.
}
\usage{
dmixpois(x, mu, sd, invlink = exp, GHorder = 5)
pmixpois(q, mu, sd, invlink = exp, lower.tail = TRUE, GHorder = 5)
qmixpois(p, mu, sd, invlink = exp, lower.tail = TRUE, GHorder = 5)
rmixpois(n, mu, sd, invlink = exp)
}
\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{mu}{
    Mean of the linear predictor. A single numeric value.
  }
  \item{sd}{
    Standard deviation of the linear predictor. A single numeric value.
  }
  \item{invlink}{
    Inverse link function. A function in the \R language,
    used to transform the linear predictor into the
    parameter \code{lambda} of the Poisson distribution.
  }
  \item{lower.tail}{
    Logical. If \code{TRUE} (the default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
  }
  \item{GHorder}{
    Number of quadrature points in the Gauss-Hermite quadrature approximation.
    A small positive integer.
  }
}
\details{
  These functions are analogous to
  \code{\link{dpois}}
  \code{\link{ppois}},
  \code{\link{qpois}} and
  \code{\link{rpois}}
  except that they apply to a mixture of Poisson distributions.

  In effect, the Poisson mean parameter \code{lambda} is randomised
  by setting \code{lambda = invlink(Z)} where \code{Z}
  has a Gaussian \eqn{N(\mu,\sigma^2)}{N(\mu, \sigma^2)} distribution.
  The default is \code{invlink=exp} which means that
  \code{lambda} is lognormal. Set \code{invlink=I} to assume
  that \code{lambda} is approximately Normal.

  For \code{dmixpois}, \code{pmixpois} and \code{qmixpois},
  the probability distribution is approximated using Gauss-Hermite
  quadrature. For \code{rmixpois}, the deviates are simulated
  exactly.
}
\value{
  Numeric vector:
  \code{dmixpois} gives probability masses,
  \code{ppois} gives cumulative probabilities,
  \code{qpois} gives (non-negative integer) quantiles, and
  \code{rpois} generates (non-negative integer) random deviates.
}
\seealso{
  \code{\link{dpois}},
  \code{\link{gauss.hermite}}.
}
\examples{
  dmixpois(7, 10, 1, invlink = I)
  dpois(7, 10)

  pmixpois(7, log(10), 0.2)
  ppois(7, 10)

  qmixpois(0.95, log(10), 0.2)
  qpois(0.95, 10)

  x <- rmixpois(100, log(10), log(1.2))
  mean(x)
  var(x)
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege
}
\keyword{distribution}
