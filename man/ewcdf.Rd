\name{ewcdf}
\alias{ewcdf}
\title{Weighted Empirical Cumulative Distribution Function}
\description{
  Compute a weighted version of the
  empirical cumulative distribution function.
}
\usage{
ewcdf(x, weights = NULL, normalise=TRUE, adjust=1)
}
\arguments{
  \item{x}{Numeric vector of observations.}
  \item{weights}{
    Optional. Numeric vector of non-negative weights for \code{x}.
    Defaults to equal weight 1 for each entry of \code{x}.
  }
  \item{normalise}{
    Logical value indicating whether the weights should be rescaled
    so that they sum to 1.
  }
  \item{adjust}{
    Numeric value. Adjustment factor.
    The weights will be multiplied by \code{adjust}.
  }
}
\details{
  This is a modification of the standard function \code{\link{ecdf}}
  allowing the observations \code{x} to have weights.

  The weighted e.c.d.f. (empirical cumulative distribution function)
  \code{Fn} is defined so that, for any real number \code{y}, the value of
  \code{Fn(y)} is equal to the total weight of all entries of
  \code{x} that are less than or equal to \code{y}. That is
  \code{Fn(y) = sum(weights[x <= y])}.

  Thus \code{Fn} is a step function which jumps at the
  values of \code{x}. The height of the jump at a point \code{y}
  is the total weight of all entries in \code{x} 
  number of tied observations at that value.  Missing values are
  ignored.

  If \code{weights} is omitted, the default is equivalent to
  \code{ecdf(x)} except for the class membership.

  The result of \code{ewcdf} is a function, of class \code{"ewcdf"},
  inheriting from the classes \code{"ecdf"} (if \code{normalise=TRUE})
  and \code{"stepfun"}.
  The class \code{ewcdf} has methods for \code{print} and \code{quantile}.
  The inherited classes \code{ecdf} and \code{stepfun}
  have methods for \code{plot} and \code{summary}.
}
\value{
  A function, of class \code{"ewcdf"}, inheriting from 
  \code{"ecdf"} (if \code{normalise=TRUE}) and \code{"stepfun"}. 
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{ecdf}}.
  
  \code{\link{quantile.ewcdf}}
}
\examples{
   x <- rnorm(100)
   w <- runif(100)
   plot(e <- ewcdf(x,w))
   e
}
\keyword{nonparametric}
\keyword{univar}
