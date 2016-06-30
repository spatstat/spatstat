\name{logLik.ppm}
\alias{logLik.ppm}
\alias{deviance.ppm}
\alias{AIC.ppm}
\alias{extractAIC.ppm}
\alias{nobs.ppm}
\title{Log Likelihood and AIC for Point Process Model}
\description{
  Extracts the log likelihood, deviance, and AIC
  of a fitted Poisson point process
  model, or analogous quantities based on the pseudolikelihood
  or logistic likelihood for a fitted Gibbs point process model.
}
\usage{
\method{logLik}{ppm}(object, \dots, new.coef=NULL, warn=TRUE, absolute=FALSE)

\method{deviance}{ppm}(object, \dots)

\method{AIC}{ppm}(object, \dots, k=2, takeuchi=TRUE)

\method{extractAIC}{ppm}(fit, scale=0, k=2, \dots, takeuchi=TRUE)

\method{nobs}{ppm}(object, \dots)
}
\arguments{
  \item{object,fit}{Fitted point process model.
    An object of class \code{"ppm"}.
  }
  \item{\dots}{Ignored.}
  \item{warn}{
    If \code{TRUE}, a warning is given when the
    pseudolikelihood or logistic likelihood
    is returned instead of the likelihood.
  }
  \item{absolute}{
    Logical value indicating whether to include 
    constant terms in the loglikelihood. 
  }
  \item{scale}{Ignored.}
  \item{k}{Numeric value specifying the weight of the
    equivalent degrees of freedom in the AIC. See Details.
  }
  \item{new.coef}{
    New values for the canonical parameters of the model.
    A numeric vector of the same length as \code{coef(object)}.
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
  \code{\link[stats]{deviance}},
  \code{\link[stats]{extractAIC}}  and
  \code{\link[stats]{nobs}}
  for the class \code{"ppm"}.

  An object of class \code{"ppm"} represents a fitted
  Poisson or Gibbs point process model.
  It is obtained from the model-fitting function \code{\link{ppm}}.
  
  The method \code{logLik.ppm} computes the 
  maximised value of the log likelihood for the fitted model \code{object}
  (as approximated by quadrature using the Berman-Turner approximation)
  is extracted. If \code{object} is not a Poisson process, the maximised log
  \emph{pseudolikelihood} is returned, with a warning (if \code{warn=TRUE}).

  The Akaike Information Criterion AIC for a fitted model is defined as
  \deqn{
    AIC = -2 \log(L) + k \times \mbox{penalty}
  }{
    AIC = -2 * log(L) + k * penalty
  }
  where \eqn{L} is the maximised likelihood of the fitted model,
  and \eqn{\mbox{penalty}}{penalty} is a penalty for model complexity,
  usually equal to the effective degrees of freedom of the model.
  The method \code{extractAIC.ppm} returns the \emph{analogous} quantity
  \eqn{AIC*} in which \eqn{L} is replaced by \eqn{L*},
  the quadrature approximation
  to the likelihood (if \code{fit} is a Poisson model)
  or the pseudolikelihood or logistic likelihood
  (if \code{fit} is a Gibbs model).

  The \eqn{\mbox{penalty}}{penalty} term is calculated
  as follows. If \code{takeuchi=FALSE} then \eqn{\mbox{penalty}}{penalty} is
  the number of fitted parameters. If \code{takeuchi=TRUE} then
  \eqn{\mbox{penalty} = \mbox{trace}(J H^{-1})}{penalty = trace(J H^(-1))}
  where \eqn{J} and \eqn{H} are the estimated variance and hessian,
  respectively, of the composite score.
  These two choices are equivalent for a Poisson process.
  
  The method \code{nobs.ppm} returns the number of points
  in the original data point pattern to which the model was fitted.
  
  The \R function \code{\link[stats]{step}} uses these methods.
}
\value{
  \code{logLik} returns a numerical value, belonging to the class
  \code{"logLik"}, with an attribute \code{"df"} giving the degrees of
  freedom.
  
  \code{AIC} returns a numerical value.

  \code{extractAIC} returns a numeric vector of length 2
  containing the degrees of freedom and the AIC value.

  \code{nobs} returns an integer value.
}
\references{
  Varin, C. and Vidoni, P. (2005)
  A note on composite likelihood inference and model selection.
  \emph{Biometrika} \bold{92}, 519--528.
}
\seealso{
  \code{\link{ppm}},  
  \code{\link{as.owin}},
  \code{\link{coef.ppm}},
  \code{\link{fitted.ppm}},
  \code{\link{formula.ppm}},
  \code{\link{model.frame.ppm}},
  \code{\link{model.matrix.ppm}},
  \code{\link{plot.ppm}},
  \code{\link{predict.ppm}},
  \code{\link{residuals.ppm}},
  \code{\link{simulate.ppm}},
  \code{\link{summary.ppm}},
  \code{\link{terms.ppm}},
  \code{\link{update.ppm}},
  \code{\link{vcov.ppm}}.
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\examples{
  data(cells)
  fit <- ppm(cells, ~x)
  nobs(fit)
  logLik(fit)
  deviance(fit)
  extractAIC(fit)
  AIC(fit)
  step(fit)
}
\keyword{spatial}
\keyword{models}

