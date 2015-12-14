\name{coef.ppm}
\alias{coef.ppm}
\title{
  Coefficients of Fitted Point Process Model
}
\description{
  Given a point process model fitted to a point pattern,
  extract the coefficients of the fitted model.
  A method for \code{coef}.
}
\usage{
  \method{coef}{ppm}(object, \dots)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"})
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
  (object of class \code{"ppm"}). Such objects are produced by the maximum
  pseudolikelihood fitting algorithm \code{\link{ppm}}).

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

  Use \code{\link{print.ppm}} to print a more useful
  description of the fitted model.
}
\seealso{
 \code{\link{print.ppm}},
 \code{\link{ppm.object}},
 \code{\link{ppm}}
}
\examples{
    data(cells)

    poi <- ppm(cells, ~1, Poisson())
    coef(poi)
    # This is the log of the fitted intensity of the Poisson process

    stra <- ppm(cells, ~1, Strauss(r=0.07))
    coef(stra)

    # The two entries "(Intercept)" and "Interaction"
    # are respectively log(beta) and log(gamma)
    # in the usual notation for Strauss(beta, gamma, r)

}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
