\name{vcov.kppm}
\alias{vcov.kppm}
\title{Variance-Covariance Matrix for a Fitted Cluster Point Process Model}
\description{
  Returns the variance-covariance matrix of the estimates of the
  parameters of a fitted cluster point process model.
}
\usage{
   \method{vcov}{kppm}(object, ...,
          what=c("vcov", "corr", "fisher", "internals"))
}
\arguments{
  \item{object}{
    A fitted cluster point process model (an object of class
    \code{"kppm"}.)
  }
  \item{\dots}{
    Ignored.
  }
  \item{what}{
    Character string (partially-matched)
    that specifies what matrix is returned.
    Options are \code{"vcov"} for the variance-covariance matrix,
    \code{"corr"} for the correlation matrix, and
    \code{"fisher"} for the Fisher information matrix.
  }
}
\details{
  This function computes the asymptotic variance-covariance
  matrix of the estimates of the canonical (regression) parameters in the
  cluster point process model \code{object}. It is a method for the 
  generic function \code{\link{vcov}}.
  
  The result is an \code{n * n} matrix where \code{n =
    length(coef(model))}.

  To calculate a confidence interval for a regression parameter,
  use \code{\link[stats]{confint}} as shown in the examples.
}
\value{
  A square matrix. 
}
\references{
  Waagepetersen, R. (2007)
  Estimating functions for inhomogeneous spatial point processes
  with incomplete covariate data.
  \emph{Biometrika} \bold{95}, 351--363.
}
\author{
  Abdollah Jalilian and Rasmus Waagepetersen.
  Ported to \pkg{spatstat} by Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
}
\seealso{
  \code{\link{kppm}}, 
  \code{\link{vcov}}, 
  \code{\link{vcov.ppm}}
}
\examples{
   data(redwood)
   fit <- kppm(redwood, ~ x + y)
   vcov(fit)
   vcov(fit, what="corr")

   # confidence interval
   confint(fit)
   # cross-check the confidence interval by hand:
   sd <- sqrt(diag(vcov(fit)))
   t(coef(fit) + 1.96 * outer(sd, c(lower=-1, upper=1)))
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
