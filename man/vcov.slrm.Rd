\name{vcov.slrm}
\alias{vcov.slrm}
\title{Variance-Covariance Matrix for a Fitted Spatial Logistic Regression}
\description{
  Returns the variance-covariance matrix of the estimates of the
  parameters of a point process model that was fitted by
  spatial logistic regression.
}
\usage{
  \method{vcov}{slrm}(object, \dots,
         what=c("vcov", "corr", "fisher", "Fisher")) 
}
\arguments{
  \item{object}{A fitted point process model of class \code{"slrm"}.}
  \item{\dots}{Ignored.}
  \item{what}{Character string (partially-matched)
    that specifies what matrix is returned.
    Options are \code{"vcov"} for the variance-covariance matrix,
    \code{"corr"} for the correlation matrix, and
    \code{"fisher"} or \code{"Fisher"} for the Fisher information matrix.
  }
}
\details{
  This function computes the asymptotic variance-covariance
  matrix of the estimates of the canonical parameters in the
  point process model \code{object}. It is a method for the 
  generic function \code{\link{vcov}}.

  \code{object} should be an object of class \code{"slrm"}, typically
  produced by \code{\link{slrm}}. It represents a Poisson point process
  model fitted by spatial logistic regression.

  The canonical parameters of the fitted model \code{object}
  are the quantities returned by \code{coef.slrm(object)}.
  The function \code{vcov} calculates the variance-covariance matrix
  for these parameters.
  
  The argument \code{what} provides three options:
  \describe{
    \item{\code{what="vcov"}}{
      return the variance-covariance matrix of the parameter estimates
    }
    \item{\code{what="corr"}}{
      return the correlation matrix of the parameter estimates
    }
    \item{\code{what="fisher"}}{
      return the observed Fisher information matrix.
    }
  }
  In all three cases, the result is a square matrix.
  The rows and columns of the matrix correspond to the canonical
  parameters given by \code{\link{coef.slrm}(object)}. The row and column
  names of the matrix are also identical to the names in
  \code{\link{coef.slrm}(object)}.

  Note that standard errors and 95\% confidence intervals for
  the coefficients can also be obtained using
  \code{confint(object)} or \code{coef(summary(object))}.

  Standard errors for the fitted intensity can be obtained
  using \code{\link{predict.slrm}}.
}
\section{Error messages}{
  An error message that reports
  \emph{system is computationally singular} indicates that the
  determinant of the Fisher information matrix was either too large 
  or too small for reliable numerical calculation.
  This can occur because of numerical overflow or
  collinearity in the covariates. 
}
\value{
  A square matrix.
}
\examples{
  X <- rpoispp(42)
  fit <- slrm(X ~ x + y)
  vcov(fit)
  vcov(fit, what="corr")
  vcov(fit, what="f")
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner \email{r.turner@auckland.ac.nz}.
}
\seealso{
  \code{\link{vcov}} for the generic,

  \code{\link{slrm}} for information about fitted models,

  \code{\link{predict.slrm}} for other kinds of calculation about the model,

  \code{\link[stats]{confint}} for confidence intervals.
}
\references{
  Baddeley, A., Berman, M., Fisher, N.I., Hardegen, A., Milne, R.K.,
  Schuhmacher, D., Shah, R. and Turner, R. (2010)
  Spatial logistic regression and change-of-support
  for spatial Poisson point processes.
  \emph{Electronic Journal of Statistics}
  \bold{4}, 1151--1201.
  {doi: 10.1214/10-EJS581}
}
\keyword{spatial}
\keyword{methods}
\keyword{models}

