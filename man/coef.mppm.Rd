\name{coef.mppm}
\alias{coef.mppm}
\title{
  Coefficients of Point Process Model Fitted to Multiple Point Patterns
}
\description{
  Given a point process model fitted to a list of point patterns,
  extract the coefficients of the fitted model.
  A method for \code{coef}.
}
\usage{
  \method{coef}{mppm}(object, \dots)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"mppm"})
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  A vector containing the fitted coefficients.
}
\details{
  This function is a method for the generic function \code{\link{coef}}.
  
  The argument \code{object} must be a fitted point process model
  (object of class \code{"mppm"}) produced by the 
  fitting algorithm \code{\link{mppm}}). This represents a
  point process model that has been fitted
  to a list of several point pattern datasets. See \code{\link{mppm}}
  for information.

  This function extracts the vector of coefficients of the fitted model.
  This is the estimate of the parameter vector
  \eqn{\theta}{theta} such that the conditional intensity of the model
  is of the form
  \deqn{
    \lambda(u,x) = \exp(\theta S(u,x))
  }{
    lambda(u,x) = exp(theta . S(u,x))
  }
  where \eqn{S(u,x)} is a (vector-valued) statistic.

  For example, if the model \code{object} is the uniform Poisson process,
  then \code{coef(object)} will yield a single value
  (named \code{"(Intercept)"}) which is the logarithm of the
  fitted intensity of the Poisson process.

  Use \code{\link{print.mppm}} to print a more useful
  description of the fitted model.
}
\seealso{
 \code{\link{print.mppm}},
 \code{\link{mppm}}
}
\examples{
    data(waterstriders)
    H <- hyperframe(X=waterstriders)
    fit.Poisson <- mppm(X ~ 1, H)
    coef(fit.Poisson)

    # The single entry "(Intercept)" 
    # is the log of the fitted intensity of the Poisson process

    fit.Strauss <- mppm(X~1, H, Strauss(7))
    coef(fit.Strauss)

    # The two entries "(Intercept)" and "Interaction"
    # are respectively log(beta) and log(gamma)
    # in the usual notation for Strauss(beta, gamma, r)

}
\author{Adrian Baddeley
  \email{adrian.baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
