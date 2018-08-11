\name{weighted.median}
\alias{weighted.median}
\alias{weighted.quantile}
\alias{weighted.var}
\title{
  Weighted Median, Quantiles or Variance
}
\description{
  Compute the median, quantiles or variance of a set of numbers which have
  weights associated with them.
}
\usage{
weighted.median(x, w, na.rm = TRUE)

weighted.quantile(x, w, probs=seq(0,1,0.25), na.rm = TRUE)

weighted.var(x, w, na.rm = TRUE)
}
\arguments{
  \item{x}{
    Data values.
    A vector of numeric values, for which the median or quantiles are required.
  }
  \item{w}{
    Weights.
    A vector of nonnegative numbers, of the same length as \code{x}.
  }
  \item{probs}{
    Probabilities for which the quantiles should be computed.
    A numeric vector of values between 0 and 1.
  }
  \item{na.rm}{
    Logical. Whether to ignore \code{NA} values.
  }
}
\details{
  The \code{i}th observation \code{x[i]} is treated as having
  a weight proportional to \code{w[i]}.

  The weighted median is a value \code{m}
  such that the total weight of data to the left of \code{m}
  is equal to half the total weight.
  If there is no such value, linear interpolation is performed.
}
\value{
  A numeric value or vector.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[stats]{quantile}}, \code{\link[stats]{median}}.
}
\examples{
  x <- 1:20
  w <- runif(20)
  weighted.median(x, w)
  weighted.quantile(x, w)
  weighted.var(x, w)
}
\keyword{math}
