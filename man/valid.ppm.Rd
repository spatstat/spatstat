\name{valid.ppm}
\alias{valid.ppm}
\title{
  Check Whether Point Process Model is Valid
}
\description{
  Determines whether a fitted point process model 
  satisfies the integrability conditions for existence of the point process.
}
\usage{
valid.ppm(object)
}
\arguments{
  \item{object}{
    Fitted point process model (object of class \code{"ppm"}).
  }
}
\details{
  The model-fitting function \code{\link{ppm}}
  fits Gibbs point process models to point pattern data.
  By default, \code{\link{ppm}} does not check whether the
  fitted model actually exists as a point process. This checking
  is done by \code{valid.ppm}.

  Unlike a regression model, which is well-defined for any values
  of the fitted regression coefficients, a Gibbs point process model
  is only well-defined if the fitted interaction parameters 
  satisfy some constraints. 
  A famous example is the Strauss process (see \code{\link{Strauss}})
  which exists only when the interaction parameter \eqn{\gamma}{gamma}
  is less than or equal to 1. For values \eqn{\gamma > 1}{gamma > 1},
  the probability density is not integrable and the process does not
  exist (and cannot be simulated).

  By default, \code{\link{ppm}} does not enforce the constraint that
  a fitted Strauss process (for example) must satisfy
  \eqn{\gamma \le 1}{gamma <= 1}.
  This is because a fitted parameter value of \eqn{\gamma > 1}{gamma > 1} 
  could be useful information for data analysis, as it indicates that
  the Strauss model is not appropriate, and suggests a clustered model should be
  fitted.

  The function \code{valid.ppm} checks whether the fitted model
  \code{object} specifies a well-defined point process. It returns
  \code{TRUE} if the model is well-defined.

  Another possible reason for invalid models is that the data may not
  be adequate for estimation of the model parameters. In this case,
  some of the fitted coefficients could be \code{NA} or infinite
  values. If this happens
  then \code{valid.ppm} returns \code{FALSE}.

  Use the function \code{\link{project.ppm}} to force the fitted model
  to be valid.
}
\value{
  Logical.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{ppm}},
  \code{\link{project.ppm}}
}
\examples{
   fit1 <- ppm(cells, ~1, Strauss(0.1))
   valid.ppm(fit1)
   fit2 <- ppm(redwood, ~1, Strauss(0.1))
   valid.ppm(fit2)
}
\keyword{spatial}
\keyword{models}
