\name{gauss.hermite}
\alias{gauss.hermite}
\title{
  Gauss-Hermite Quadrature Approximation to Expectation for Normal Distribution
}
\description{
  Calculates an approximation to the expected value of any function of a
  normally-distributed random variable, using Gauss-Hermite quadrature.
}
\usage{
gauss.hermite(f, mu = 0, sd = 1, ..., order = 5)
}
\arguments{
  \item{f}{
    The function whose moment should be approximated.
  }
  \item{mu}{
    Mean of the normal distribution.
  }
  \item{sd}{
    Standard deviation of the normal distribution.
  }
  \item{\dots}{
    Additional arguments passed to \code{f}.
  }
  \item{order}{
    Number of quadrature points in the Gauss-Hermite quadrature
    approximation. A small positive integer.
  }
}
\details{
  This algorithm calculates the approximate expected value of
  \code{f(Z)} when \code{Z} is a normally-distributed random
  variable with mean \code{mu} and standard deviation \code{sd}.
  The expected value is an integral with respect to the
  Gaussian density; this integral is approximated
  using Gauss-Hermite quadrature.

  The argument \code{f} should be a function in the \R language
  whose first argument is the variable \code{Z}. Additional arguments
  may be passed through \code{\dots}. The value returned by \code{f}
  may be a single numeric value, a vector, or a matrix. The values
  returned by \code{f} for different values of \code{Z} must have
  compatible dimensions.

  The result is a weighted average of several values of \code{f}.
}
\value{
  Numeric value, vector or matrix.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  ,
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}.
}
\examples{
  gauss.hermite(function(x) x^2, 3, 1)
}
\keyword{math}
