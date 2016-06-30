\name{quantile.density}
\alias{quantile.density}
\title{
  Quantiles of a Density Estimate
}
\description{
  Given a kernel estimate of a probability density, compute quantiles.
}
\usage{
\method{quantile}{density}(x, probs = seq(0, 1, 0.25), names = TRUE,
        \dots, warn = TRUE)
}
\arguments{
  \item{x}{
    Object of class \code{"density"} computed by a method for
    \code{\link[stats]{density}}
  }
  \item{probs}{
    Numeric vector of probabilities for which the quantiles are required.
  }
  \item{names}{
    Logical value indicating whether to attach names (based on
    \code{probs}) to the result.
  }
  \item{\dots}{
    Ignored.
  }
  \item{warn}{
    Logical value indicating whether to issue a warning if the
    density estimate \code{x} had to be renormalised because it
    was computed in a restricted interval.
  }
}
\details{
  This function calculates quantiles of the probability distribution
  whose probability density has been estimated and stored in the object
  \code{x}. The object \code{x} must belong to the class \code{"density"},
  and would typically have been obtained from a call to the function
  \code{\link[stats]{density}}.

  The probability density is first normalised so that the total
  probability is equal to 1. A warning is issued if the density
  estimate was restricted to an interval (i.e. if \code{x}
  was created by a call to   \code{\link[stats]{density}} which
  included either of the arguments \code{from} and \code{to}).

  Next, the density estimate is numerically integrated to obtain an estimate
  of the cumulative distribution function \eqn{F(x)}. Then
  for each desired probability \eqn{p}, the algorithm finds the
  corresponding quantile \eqn{q}.

  The quantile \eqn{q} corresponding to probability \eqn{p}
  satisfies \eqn{F(q) = p} up to
  the resolution of the grid of values contained in \code{x}.
  The quantile is computed from the right, that is,
  \eqn{q} is the smallest available value of \eqn{x} such that
  \eqn{F(x) \ge p}{F(x) >= p}.
}
\value{
  A numeric vector containing the quantiles.
}
\author{
  \adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link[stats]{quantile}}, 
  \code{\link{quantile.ewcdf}}, 
  \code{\link{quantile.im}},
  \code{\link{CDF}}.
}
\examples{
   dd <- density(runif(10))
   quantile(dd)
}
\keyword{methods}
\keyword{univar}
\keyword{nonparametric}
