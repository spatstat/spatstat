\name{logLik.mppm}
\alias{logLik.mppm}
\alias{AIC.mppm}
\alias{extractAIC.mppm}
\alias{nobs.mppm}
\alias{getCall.mppm}
\alias{terms.mppm}
\title{Log Likelihood and AIC for Multiple Point Process Model}
\description{
  For a point process model that has been fitted to multiple point
  patterns, these functions extract the log likelihood and AIC,
  or analogous quantities based on the pseudolikelihood.
}
\usage{
 \method{logLik}{mppm}(object, \dots, warn=TRUE)

 \method{AIC}{mppm}(object, \dots, k=2, takeuchi=TRUE)

 \method{extractAIC}{mppm}(fit, scale = 0, k = 2, \dots, takeuchi = TRUE) 

 \method{nobs}{mppm}(object, \dots)

 \method{getCall}{mppm}(x, \dots)

 \method{terms}{mppm}(x, \dots)
}
\arguments{
  \item{object,fit,x}{
    Fitted point process model (fitted to multiple point
    patterns). An object of class \code{"mppm"}.
  }
  \item{\dots}{Ignored.}
  \item{warn}{
    If \code{TRUE}, a warning is given when the
    pseudolikelihood is returned instead of the likelihood.
  }
  \item{scale}{Ignored.}
  \item{k}{Numeric value specifying the weight of the
    equivalent degrees of freedom in the AIC. See Details.
  }
  \item{takeuchi}{
    Logical value specifying whether to use the Takeuchi penalty
    (\code{takeuchi=TRUE}) or the
    number of fitted parameters (\code{takeuchi=FALSE})
    in calculating AIC.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link[stats]{logLik}},
  \code{\link[stats]{AIC}},
  \code{\link[stats]{extractAIC}},
  \code{\link[stats]{terms}} and
  \code{\link[stats:update]{getCall}} 
  for the class \code{"mppm"}.

  An object of class \code{"mppm"} represents a fitted
  Poisson or Gibbs point process model fitted to several point patterns.
  It is obtained from the model-fitting function \code{\link{mppm}}.

  The method \code{logLik.mppm} extracts the 
  maximised value of the log likelihood for the fitted model
  (as approximated by quadrature using the Berman-Turner approximation).
  If \code{object} is not a Poisson process, the maximised log
  \emph{pseudolikelihood} is returned, with a warning.

  The Akaike Information Criterion AIC for a fitted model is defined as
  \deqn{
    AIC = -2 \log(L) + k \times \mbox{penalty}
  }{
    AIC = -2 * log(L) + k * penalty
  }
  where \eqn{L} is the maximised likelihood of the fitted model,
  and \eqn{\mbox{penalty}}{penalty} is a penalty for model complexity,
  usually equal to the effective degrees of freedom of the model.
  The method \code{extractAIC.mppm} returns the \emph{analogous} quantity
  \eqn{AIC*} in which \eqn{L} is replaced by \eqn{L*},
  the quadrature approximation
  to the likelihood (if \code{fit} is a Poisson model)
  or the pseudolikelihood (if \code{fit} is a Gibbs model).

  The \eqn{\mbox{penalty}}{penalty} term is calculated
  as follows. If \code{takeuchi=FALSE} then \eqn{\mbox{penalty}}{penalty} is
  the number of fitted parameters. If \code{takeuchi=TRUE} then
  \eqn{\mbox{penalty} = \mbox{trace}(J H^{-1})}{penalty = trace(J H^(-1))}
  where \eqn{J} and \eqn{H} are the estimated variance and hessian,
  respectively, of the composite score.
  These two choices are equivalent for a Poisson process.
  
  The method \code{nobs.mppm} returns the total number of points
  in the original data point patterns to which the model was fitted.

  The method \code{getCall.mppm} extracts the original call to
  \code{\link{mppm}} which caused the model to be fitted.
  
  The method \code{terms.mppm} extracts the covariate terms in the
  model formula as a \code{terms} object. Note that these terms do not
  include the interaction component of the model.
  
  The \R function \code{\link[stats]{step}} uses these methods.
}
\value{
  See the help files for the corresponding generic functions.
}
\seealso{
  \code{\link{mppm}}
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\author{
  Adrian Baddeley, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  \spatstatAuthors.
}
\examples{
  fit <- mppm(Bugs ~ x, hyperframe(Bugs=waterstriders))
  logLik(fit)
  AIC(fit)
  nobs(fit)
  getCall(fit)
}
\keyword{spatial}
\keyword{models}

