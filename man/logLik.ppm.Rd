\name{logLik.ppm}
\alias{logLik.ppm}
\alias{extractAIC.ppm}
\alias{nobs.ppm}
\title{Log Likelihood and AIC for Point Process Model}
\description{
  Extracts the log likelihood, deviance, and AIC
  of a fitted Poisson point process
  model, or analogous quantities based on the pseudolikelihood for
  a fitted Gibbs point process model.
}
\usage{
\method{logLik}{ppm}(object, ..., warn=TRUE)
\method{extractAIC}{ppm}(fit, scale=0, k=2, \dots)
\method{nobs}{ppm}(object, ...)
}
\arguments{
  \item{object,fit}{Fitted point process model.
    An object of class \code{"ppm"}.
  }
  \item{\dots}{Ignored.}
  \item{warn}{
    If \code{TRUE}, a warning is given when the
    pseudolikelihood is returned instead of the likelihood.
  }
  \item{scale}{Ignored.}
  \item{k}{Numeric value specifying the weight of the
    equivalent degrees of freedom in the AIC. See Details.}
}
\details{
  These functions are methods for the generic commands
  \code{\link{logLik}},
  \code{\link{extractAIC}}  and
  \code{\link{nobs}}
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
    AIC = -2 \log(L) + k \times \mbox{edf}
  }{
    AIC = -2 * log(L) + k * edf
  }
  where \eqn{L} is the maximised likelihood of the fitted model,
  and \eqn{\mbox{edf}}{edf} is the effective degrees of freedom
  of the model.
  The method \code{extractAIC.ppm} returns the \emph{analogous} quantity
  \eqn{AIC*} in which \eqn{L} is replaced by \eqn{L*},
  the quadrature approximation
  to the likelihood (if \code{fit} is a Poisson model)
  or the pseudolikelihood (if \code{fit} is a Gibbs model).

  The method \code{nobs.ppm} returns the number of points
  in the original data point pattern to which the model was fitted.
  
  The \R functions \code{\link{AIC}} and \code{\link{step}} use
  these methods.
}
\value{
  A numerical value.
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
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\examples{
  data(cells)
  fit <- ppm(cells, ~x)
  nobs(fit)
  logLik(fit)
  extractAIC(fit)
  AIC(fit)
  step(fit)
}
\keyword{spatial}
\keyword{models}

