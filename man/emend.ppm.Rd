\name{emend.ppm}
\alias{emend.ppm}
\alias{project.ppm}
\title{
  Force Point Process Model to be Valid
}
\description{
  Ensures that a fitted point process model 
  satisfies the integrability conditions for existence of the point process.
}
\usage{
project.ppm(object, \dots, fatal=FALSE, trace=FALSE)

\method{emend}{ppm}(object, \dots, fatal=FALSE, trace=FALSE)
}
\arguments{
  \item{object}{
    Fitted point process model (object of class \code{"ppm"}).
  }
  \item{\dots}{Ignored.}
  \item{fatal}{
    Logical value indicating whether to generate an error
    if the model cannot be projected to a valid model.
  }
  \item{trace}{
    Logical value indicating whether to print a trace
    of the decision process. 
  }
}
\details{
  The functions \code{emend.ppm} and \code{project.ppm} are identical:
  \code{emend.ppm} is a method for the generic \code{\link{emend}},
  while \code{project.ppm} is an older name for the same function.

  The purpose of the function is to ensure that a fitted model
  is valid.
  
  The model-fitting function \code{\link{ppm}}
  fits Gibbs point process models to point pattern data.
  By default, the fitted model returned by \code{\link{ppm}} may not
  actually exist as a point process. 

  First, some of the fitted coefficients of the model
  may be \code{NA} or infinite values. 
  This usually occurs when the data are insufficient to estimate
  all the parameters. The model is said to be
  \emph{unidentifiable} or \emph{confounded}.

  Second, unlike a regression model, which is well-defined for any finite values
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

  The function \code{emend.ppm} or \code{project.ppm}
  modifies the model \code{object}
  so that the model is valid. It 
  identifies the terms in the model \code{object}
  that are associated with illegal parameter values (i.e. parameter
  values which are either \code{NA}, infinite, or outside their permitted
  range). It considers all possible sub-models of \code{object}
  obtained by deleting one or more
  of these terms. It identifies which of these submodels are valid,
  and chooses the valid submodel with the largest pseudolikelihood. The result
  of \code{emend.ppm} or \code{project.ppm} is the
  true maximum pseudolikelihood fit to the data.

  For large datasets or complex models, the algorithm used in
  \code{emend.ppm} or 
  \code{project.ppm} may be time-consuming, because it takes time to
  compute all the sub-models. A faster, approximate
  algorithm can be applied by setting
  \code{spatstat.options(project.fast=TRUE)}. This produces a
  valid submodel, which may not be the maximum pseudolikelihood submodel.

  Use the function \code{\link{valid.ppm}} to check whether a fitted model
  object specifies a well-defined point process.

  Use the expression \code{all(is.finite(coef(object)))} to determine
  whether all parameters are identifiable.
}
\value{
  Another point process model (object of class \code{"ppm"}).
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{ppm}},
  \code{\link{valid.ppm}},
  \code{\link{emend}},
  \code{\link{spatstat.options}}
}
\examples{
   fit <- ppm(redwood, ~1, Strauss(0.1))
   coef(fit)
   fit2 <- emend(fit)
   coef(fit2)
}
\keyword{spatial}
\keyword{models}
