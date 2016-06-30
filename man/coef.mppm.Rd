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
  Either a vector containing the fitted coefficients,
  or a data frame containing the fitted coefficients for each point pattern.
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

  If the fitted model includes random effects (i.e. if the argument
  \code{random} was specified in the call to \code{\link{mppm}}),
  then the fitted coefficients are different for each point pattern
  in the original data, so \code{coef(object)} is a data frame
  with one row for each point pattern, and one column for each
  parameter. Use \code{\link{fixef.mppm}} to extract the vector of fixed effect
  coefficients, and \code{\link{ranef.mppm}} to extract the random
  effect coefficients at each level.

  Use \code{\link{print.mppm}} to print a more useful
  description of the fitted model.
}
\seealso{
  \code{\link{fixef.mppm}} and \code{\link{ranef.mppm}}
  for the fixed and random effect coefficients in a model that includes
  random effects.
  
 \code{\link{print.mppm}},
 \code{\link{mppm}}
}
\examples{
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

    # Tweak data to exaggerate differences
    H$X[[1]] <- rthin(H$X[[1]], 0.3)
    # Model with random effects
    fitran <- mppm(X ~ 1, H, random=~1|id)
    coef(fitran)
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\author{
  Adrian Baddeley, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented in \pkg{spatstat} by
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
