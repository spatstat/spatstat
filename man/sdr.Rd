\name{sdr}
\alias{sdr}
\title{
  Sufficient Dimension Reduction
}
\description{
  Given a point pattern and a set of predictors,
  find a minimal set of new predictors, each constructed as
  a linear combination of the original predictors.
}
\usage{
sdr(X, covariates, method = c("DR", "NNIR", "SAVE", "SIR", "TSE"),
    Dim1 = 1, Dim2 = 1, predict=FALSE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{covariates}{
    A list of pixel images (objects of class \code{"im"})
    to serve as predictor variables.
  }
  \item{method}{
    Character string indicating which method to use. See Details.
  }
  \item{Dim1}{
    Dimension of the first order Central Intensity Subspace
    (applicable when \code{method} is \code{"DR"}, \code{"NNIR"},
    \code{"SAVE"} or \code{"TSE"}). 
  }
  \item{Dim2}{
    Dimension of the second order Central Intensity Subspace
    (applicable when \code{method="TSE"}).
  }
  \item{predict}{
    Logical value indicating whether to compute the new predictors
    as well.
  }
}
\details{
  Given a point pattern \eqn{X} and predictor variables
  \eqn{Z_1, \dots, Z_p}{Z[1], ..., Z[p]},
  Sufficient Dimension Reduction methods
  (Guan and Wang, 2010) attempt to find a minimal set
  of new predictor variables, each constructed by taking a linear combination
  of the original predictors, which explain the dependence of
  \eqn{X} on   \eqn{Z_1, \dots, Z_p}{Z[1], ..., Z[p]}.
  The methods do not assume any particular form of dependence
  of the point pattern on the predictors.
  The predictors are assumed to
  be Gaussian random fields.

  Available methods are:
  \tabular{ll}{
    \code{method="DR"} \tab directional regression \cr
    \code{method="NNIR"} \tab nearest neighbour inverse regression \cr
    \code{method="SAVE"} \tab sliced average variance estimation \cr
    \code{method="SIR"} \tab sliced inverse regression \cr
    \code{method="TSE"} \tab two-step estimation \cr
  }

  The result includes a matrix \code{B} whose columns are estimates
  of the basis vectors of the space of new predictors. That is,
  the \code{j}th column of \code{B} expresses the \code{j}th new
  predictor as a linear combination of the original predictors.

  If \code{predict=TRUE}, the new predictors are also evaluated.
  They can also be evaluated using \code{\link{sdrPredict}}.
}
\value{
  A list with components \code{B, M}
  or \code{B, M1, M2} where
  \code{B} is a matrix whose columns are estimates of the basis vectors
  for the space, and \code{M} or \code{M1,M2} are matrices containing
  estimates of the kernel.

  If \code{predict=TRUE}, the result also includes a component
  \code{Y} which is a list of pixel images giving the values of the
  new predictors.
}
\examples{
   A <- sdr(bei, bei.extra, predict=TRUE)
   A
   Y1 <- A$Y[[1]]
   plot(Y1)
   points(bei, pch=".", cex=2)
   # investigate likely form of dependence
   plot(rhohat(bei, Y1))
}
\seealso{
  \code{\link{sdrPredict}} to compute the new predictors from the
  coefficient matrix.
  
  \code{\link{dimhat}} to estimate the subspace dimension.

  \code{\link{subspaceDistance}}
}
\references{
  Guan, Y. and Wang, H. (2010)
  Sufficient dimension reduction for spatial point
  processes directed by Gaussian random fields.
  \emph{Journal of the Royal Statistical Society, Series B},
  \bold{72}, 367--387.
}
\author{
  Matlab original by Yongtao Guan,
  translated to \R by Suman Rakshit.
}
\keyword{spatial}
\keyword{multivariate}

