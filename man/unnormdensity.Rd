\name{unnormdensity}
\alias{unnormdensity}
\title{
  Weighted kernel smoother
}
\description{
  An unnormalised version of kernel density estimation
  where the weights are not required to sum to 1.
  The weights may be positive, negative or zero.
}
\usage{
unnormdensity(x, ..., weights = NULL)
}
\arguments{
  \item{x}{
    Numeric vector of data
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.default}}.
    Arguments must be \emph{named}.
  }`
  \item{weights}{
    Optional numeric vector of weights for the data.
  }
}
\details{
  This is an alternative to the standard \R kernel density estimation function
  \code{\link{density.default}}.

  The standard \code{\link{density.default}}
  requires the \code{weights} to be nonnegative numbers that add up to 1,
  and returns a probability density (a function that integrates to 1).
  
  This function \code{unnormdensity} does not impose any requirement
  on the \code{weights} except that they be finite. Individual weights may be
  positive, negative or zero. The result is a function that does not
  necessarily integrate to 1 and may be negative. The result is
  the convolution of the kernel \eqn{k} with the weighted data,
  \deqn{
    f(x) = \sum_i w_i k(x- x_i)
  }{
    f(x) = sum of w[i] * k(x - x[i])
  }
  where \eqn{x_i}{x[i]} are the data points and \eqn{w_i}{w[i]} are the
  weights.

  The algorithm first selects the kernel bandwidth by
  applying \code{\link{density.default}} to the data
  \code{x} with normalised, positive weight vector
  \code{w = abs(weights)/sum(abs(weights))} and
  extracting the selected bandwidth.
  Then the result is computed by applying
  applying \code{\link{density.default}} to \code{x} twice
  using the normalised positive and negative parts of the weights.
  
  Note that the arguments \code{\dots} must be passed by name,
  i.e. in the form (\code{name=value}). Arguments that do not match
  an argument of \code{\link{density.default}} will be ignored
  \emph{silently}.
}
\value{
  Object of class \code{"density"} as described in
  \code{\link{density.default}}.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
 \code{\link{density.default}}
}
\examples{
  d <- unnormdensity(1:3, weights=c(-1,0,1))
  if(interactive()) plot(d)
}
\keyword{smooth}
