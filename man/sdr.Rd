\name{sdr}
\alias{sdr}
\title{
  Sufficient Dimension Reduction
}
\description{
  Applies a chosen algorithm for Sufficient Dimension Reduction
  to \pkg{spatstat} data objects.
}
\usage{
sdr(X, covariates, method = c("DR", "NNIR", "SAVE", "SIR", "TSE"),
    Dim1 = 1, Dim2 = 1)
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
    Character string indicating which method to use.
  }
  \item{Dim1}{
    Dimension of the Central Subspace for the first order intensity
    (applicable when \code{method} is \code{"DR"}, \code{"NNIR"},
    \code{"SAVE"} or \code{"TSE"}). 
  }
  \item{Dim2}{
    Dimension of the Central Subspace for the second moments
    (applicable when \code{method="TSE"}).
  }
}
\details{
  Insert an explanation here!
}
\value{
  A list with components \code{B, M}
  or \code{B, M1, M2} where
  \code{B} is a matrix whose columns are estimates of the basis vectors
  for the space, and \code{M} or \code{M1,M2} are matrices containing
  estimates of the kernel.
}
\examples{
   sdr(bei, bei.extra)
}
\seealso{
  \code{\link{dimhat}}, \code{\link{subspaceDistance}}
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

