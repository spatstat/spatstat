\name{vcov.kppm}
\alias{vcov.kppm}
\title{Variance-Covariance Matrix for a Fitted Cluster Point Process Model}
\description{
  Returns the variance-covariance matrix of the estimates of the
  parameters of a fitted cluster point process model.
}
\usage{
   \method{vcov}{kppm}(object, ...,
          what=c("vcov", "corr", "fisher", "internals"),
          fast = NULL, rmax = NULL, eps.rmax = 0.01,
          verbose = TRUE)
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
  \item{fast}{
    Logical specifying whether tapering (using sparse matrices from
    \pkg{Matrix}) should be used to speed up calculations. Warning: This
    is expected to underestimate the true asymptotic
    variances/covariances.
  }  
  \item{rmax}{
    Optional. The dependence range. Not usually specified by the
    user. Only used when \code{fast=TRUE}.
  }
  \item{eps.rmax}{
    Numeric. A small positive number which is used to determine \code{rmax}
    from the tail behaviour of the pair correlation function when fast
    option (\code{fast=TRUE}) is used. Namely
    \code{rmax} is the smallest value of \eqn{r}
    at which \eqn{(g(r)-1)/(g(0)-1)}
    falls below \code{eps.rmax}.
    Only used when \code{fast=TRUE}.
    Ignored if \code{rmax} is provided.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports
    during very long calculations.
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
  and Ege Rubak.
}
\seealso{
  \code{\link{kppm}}, 
  \code{\link{vcov}}, 
  \code{\link{vcov.ppm}}
}
\examples{
   data(redwood)
   fit <- kppm(redwood ~ x + y)
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
