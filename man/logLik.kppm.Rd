\name{logLik.kppm}
\alias{logLik.kppm}
\alias{AIC.kppm}
\alias{extractAIC.kppm}
\alias{nobs.kppm}
\title{Log Likelihood and AIC for Fitted Cox or Cluster Point Process Model}
\description{
  Extracts the log composite likelihood, deviance, and AIC
  of a fitted Cox or cluster point process model.
}
\usage{
\method{logLik}{kppm}(object, ...)
\method{AIC}{kppm}(object, \dots, k=2)
\method{extractAIC}{kppm}(fit, scale=0, k=2, \dots)
\method{nobs}{kppm}(object, ...)
}
\arguments{
  \item{object,fit}{Fitted point process model.
    An object of class \code{"kppm"}.
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
  for the class \code{"kppm"}.

  An object of class \code{"kppm"} represents a fitted
  Cox or cluster point process model.
  It is obtained from the model-fitting function \code{\link{kppm}}.

  These methods apply only when the model was fitted
  by maximising a composite likelihood:
  either the Palm likelihood (Tanaka et al, 2008)
  or the second order composite likelihood (Guan, 2006),
  by calling \code{\link{kppm}} with the argument \code{method="palm"}
  or \code{method="clik2"} respectively. 
  
  The method \code{logLik.kppm} computes the 
  maximised value of the log composite likelihood
  for the fitted model \code{object}.

  The methods \code{AIC.kppm} and \code{extractAIC.kppm} compute the
  Akaike Information Criterion AIC for the fitted model
  based on the composite likelihood
  \deqn{
    AIC = -2 \log(CL) + k \times \mbox{edf}
  }{
    AIC = -2 * log(CL) + k * edf
  }
  where \eqn{CL} is the maximised composite likelihood of the fitted model,
  and \eqn{\mbox{edf}}{edf} is the effective degrees of freedom
  of the model. 

  The method \code{nobs.kppm} returns the number of points
  in the original data point pattern to which the model was fitted.
  
  The \R function \code{\link{step}} uses these methods.
}
\section{Model comparison}{
  The values of log-likelihood and AIC returned by these functions
  are based on the \emph{composite likelihood} of the cluster process or Cox
  process model. They are available only when the model was fitted
  using \code{method="palm"} or \code{method="clik2"}.

  For model comparison and model selection,
  it is valid to compare the \code{logLik} values,
  or to compare the \code{AIC} values, but only when
  all the models are of class \code{"kppm"}
  and were fitted using the same \code{method}.

  For \code{method="palm"} some theoretical justification was provided by
  Tanaka et al (2008).
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
  Guan, Y. (2006) 
  A composite likelihood approach in fitting spatial point process
  models.
  \emph{Journal of the American Statistical Association}
  \bold{101}, 1502--1512.

  Tanaka, U. and Ogata, Y. and Stoyan, D. (2008) 
  Parameter estimation and model selection for
  Neyman-Scott point processes.
  \emph{Biometrical Journal} \bold{50}, 43--57.
}
\seealso{
  \code{\link{kppm}}, 
  \code{\link{logLik.ppm}}
}
\author{
  \spatstatAuthors.
}
\examples{
  fit <- kppm(redwood ~ x, "Thomas", method="palm")
  nobs(fit)
  logLik(fit)
  extractAIC(fit)
  AIC(fit)
  step(fit)
}
\keyword{spatial}
\keyword{models}

