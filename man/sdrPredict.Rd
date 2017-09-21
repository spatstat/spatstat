\name{sdrPredict}
\alias{sdrPredict}
\title{
  Compute Predictors from Sufficient Dimension Reduction
}
\description{
  Given the result of a Sufficient Dimension Reduction method,
  compute the new predictors.
}
\usage{
  sdrPredict(covariates, B)
}
\arguments{
  \item{covariates}{
    A list of pixel images (objects of class \code{"im"}).
  }
  \item{B}{
    Either a matrix of coefficients for the covariates, or the result of
    a call to \code{\link{sdr}}.
  }
}
\details{
  This function assumes that \code{\link{sdr}} has already been used to
  find a minimal set of predictors based on the \code{covariates}.
  The argument \code{B} should be either the result of \code{\link{sdr}}
  or the coefficient matrix returned as one of the
  results of \code{\link{sdr}}. The columns of this matrix define linear
  combinations of the \code{covariates}. This function evaluates those
  linear combinations, and returns a list of pixel images containing the
  new predictors.
}
\value{
  A list of pixel images  (objects of class \code{"im"})
  with one entry for each column of \code{B}.
}
\author{
  \adrian
}
\seealso{
  \code{\link{sdr}}
}
\examples{
   A <- sdr(bei, bei.extra)
   Y <- sdrPredict(bei.extra, A)
   Y
}
\keyword{spatial}
\keyword{nonparametric}
