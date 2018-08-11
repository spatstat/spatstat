\name{dimhat}
\alias{dimhat}
\title{
  Estimate Dimension of Central Subspace 
}
\description{
  Given the kernel matrix that characterises a central subspace,
  this function estimates the dimension of the subspace.
}
\usage{
  dimhat(M)
}
\arguments{
  \item{M}{
    Kernel of subspace. A symmetric, non-negative definite, numeric
    matrix, typically obtained from \code{\link{sdr}}.
  }
}
\details{
  This function computes the maximum descent estimate of
  the dimension of the central subspace with a given kernel matrix \code{M}.

  The matrix \code{M} should be the kernel matrix of a central subspace,
  which can be obtained from \code{\link{sdr}}. It must be a symmetric,
  non-negative-definite, numeric matrix.

  The algorithm finds the eigenvalues
  \eqn{\lambda_1 \ge \ldots \ge \lambda_n}{lambda[1] \ge  ...\ge lambda[n]}
  of \eqn{M},
  and then determines the index \eqn{k} for which
  \eqn{\lambda_k/\lambda_{k-1}}{lambda[k]/lambda[k-1]} is greatest.
}
\value{
  A single integer giving the estimated dimension.
}
\seealso{
  \code{\link{sdr}}, \code{\link{subspaceDistance}}
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
\keyword{array}
\keyword{algebra}
\keyword{multivariate}

