\name{logLik.dppm}
\alias{logLik.dppm}
\alias{AIC.dppm}
\alias{extractAIC.dppm}
\alias{nobs.dppm}
\title{Log Likelihood and AIC for Fitted Determinantal Point Process Model}
\description{
  Extracts the log Palm likelihood, deviance, and AIC
  of a fitted determinantal point process model.
}
\usage{
\method{logLik}{dppm}(object, ...)
\method{AIC}{dppm}(object, \dots, k=2)
\method{extractAIC}{dppm}(fit, scale=0, k=2, \dots)
\method{nobs}{dppm}(object, ...)
}
\arguments{
  \item{object,fit}{Fitted point process model.
    An object of class \code{"dppm"}.
  }
  \item{\dots}{Ignored.}
  \item{scale}{Ignored.}
  \item{k}{Numeric value specifying the weight of the
    equivalent degrees of freedom in the AIC. See Details.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link{logLik}},
  \code{\link{extractAIC}}  and
  \code{\link{nobs}}
  for the class \code{"dppm"}.

  An object of class \code{"dppm"} represents a fitted
  Cox or cluster point process model.
  It is obtained from the model-fitting function \code{\link{dppm}}.

  These methods apply only when the model was fitted
  by maximising the Palm likelihood (Tanaka et al, 2008)
  by calling \code{\link{dppm}} with the argument \code{method="palm"}.

  The method \code{logLik.dppm} computes the
  maximised value of the log Palm likelihood for the fitted model \code{object}.

  The methods \code{AIC.dppm} and \code{extractAIC.dppm} compute the
  Akaike Information Criterion AIC for the fitted model
  based on the Palm likelihood (Tanaka et al, 2008)
  \deqn{
    AIC = -2 \log(PL) + k \times \mbox{edf}
  }{
    AIC = -2 * log(PL) + k * edf
  }
  where \eqn{PL} is the maximised Palm likelihood of the fitted model,
  and \eqn{\mbox{edf}}{edf} is the effective degrees of freedom
  of the model.

  The method \code{nobs.dppm} returns the number of points
  in the original data point pattern to which the model was fitted.

  The \R function \code{\link{step}} uses these methods, but it does
  not work for determinantal models yet due to a missing implementation
  of \code{update.dppm}.
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
  Tanaka, U. and Ogata, Y. and Stoyan, D. (2008)
  Parameter estimation and model selection for
  Neyman-Scott point processes.
  \emph{Biometrical Journal} \bold{50}, 43--57.
}
\seealso{
  \code{\link{dppm}},
  \code{\link{logLik.ppm}}
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\examples{
  fit <- dppm(swedishpines ~ x, dppGauss(), method="palm")
  nobs(fit)
  logLik(fit)
  extractAIC(fit)
  AIC(fit)
}
\keyword{spatial}
\keyword{models}

