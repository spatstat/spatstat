\name{logLik.slrm}
\Rdversion{1.1}
\alias{logLik.slrm}
\title{
  Loglikelihood of Spatial Logistic Regression
}
\description{
  Computes the (maximised) loglikelihood of a fitted
  Spatial Logistic Regression model.
}
\usage{
  \method{logLik}{slrm}(object, ..., adjust = TRUE)
}
\arguments{
  \item{object}{
    a fitted spatial logistic regression model.
    An object of class \code{"slrm"}.
  }
  \item{\dots}{
    Ignored.
  }
  \item{adjust}{
    Logical value indicating whether to adjust the loglikelihood
    of the model to make it comparable with a point process
    likelihood. See Details.
  }
}
\details{
  This is a method for \code{\link[stats]{logLik}} for fitted spatial logistic
  regression models (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}). It computes the log-likelihood
  of a fitted spatial logistic regression model.

  If \code{adjust=FALSE}, the loglikelihood is computed
  using the standard formula for the loglikelihood of a
  logistic regression model for a finite set of (pixel) observations.

  If \code{adjust=TRUE} then the loglikelihood is adjusted so that it
  is approximately comparable with the likelihood of a point process
  in continuous space, by subtracting the value
  \eqn{n \log(a)}{n *  log(a)}
  where \eqn{n} is the number of points in the original point pattern
  dataset, and \eqn{a} is the area of one pixel.
}
\value{
  A numerical value.
}
\seealso{
  \code{\link{slrm}}
}
\examples{
  X <- rpoispp(42)
  fit <- slrm(X ~ x+y)
  logLik(fit)
  logLik(fit, adjust=FALSE)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  \email{adrian@maths.uwa.edu.au}
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
