\name{distcdf}
\alias{distcdf}
\title{Distribution Function of Interpoint Distance }
\description{
  Computes the cumulative distribution function of the distance
  between two independent random points in a given window
  or windows.
}
\usage{
  distcdf(W, V=W, \dots, dW=1, dV=dW, nr=1024, regularise=TRUE)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}) containing the
    first random point.
  }
  \item{V}{
    Optional. Another window containing the second random point.
    Defaults to \code{W}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the
    pixel resolution for the calculation.
  }
  \item{dV, dW}{
    Optional. Probability densities (not necessarily normalised)
    for the first and second random points respectively.
    Data in any format acceptable
    to \code{\link{as.im}}, for example, a \code{function(x,y)}
    or a pixel image or a numeric value. The default
    corresponds to a uniform distribution over the window.
  }
  \item{nr}{
    Integer. The number of values of interpoint distance \eqn{r}
    for which the CDF will be computed.
    Should be a large value!
  }
  \item{regularise}{
    Logical value indicating whether to smooth the results
    for very small distances, to avoid discretisation artefacts.
  }
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.
}
\details{
  This command computes the Cumulative Distribution Function
  \eqn{
    CDF(r) = Prob(T \le r)
  }{
    CDF(r) = Prob(T \le r)
  }
  of the Euclidean distance \eqn{T = \|X_1 - X_2\|}{T = |X1-X2|}
  between two independent random points \eqn{X_1}{X1} and \eqn{X_2}{X2}.

  In the simplest case, the command \code{distcdf(W)}, the random points are 
  assumed to be uniformly distributed in the same
  window \code{W}.

  Alternatively the two random points may be 
  uniformly distributed in two different windows \code{W} and \code{V}.

  In the most general case the first point \eqn{X_1}{X1} is random
  in the window \code{W} with a probability density proportional to
  \code{dW}, and the second point \eqn{X_2}{X2} is random in
  a different window \code{V} with probability density proportional
  to \code{dV}. The values of \code{dW} and \code{dV} must be
  finite and nonnegative.

  The calculation is performed by numerical integration of the set covariance
  function \code{\link{setcov}} for uniformly distributed points, and
  by computing the covariance function \code{\link{imcov}} in the
  general case. The accuracy of the result depends on
  the pixel resolution used to represent the windows: this is controlled
  by the arguments \code{\dots} which are passed to \code{\link{as.mask}}.
  For example use \code{eps=0.1} to specify pixels of size 0.1 units.

  The arguments \code{W} or \code{V} may also be point patterns
  (objects of class \code{"ppp"}).
  The result is the cumulative distribution function
  of the distance from a randomly selected point in the point pattern,
  to a randomly selected point in the other point pattern or window.

  If \code{regularise=TRUE} (the default), values of the cumulative
  distribution function for very short distances are smoothed to avoid
  discretisation artefacts. Smoothing is applied to all distances
  shorter than the width of 7 pixels. 
}
\seealso{
  \code{\link{setcov}},
  \code{\link{as.mask}}.
}
\examples{
 # The unit disc
 B <- disc()
 plot(distcdf(B))
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}

 
 
