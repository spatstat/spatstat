\name{matrixpower}
\alias{matrixpower}
\alias{matrixsqrt}
\alias{matrixinvsqrt}
\title{
  Power of a Matrix
}
\description{
  Evaluate a specified power of a matrix.
}
\usage{
  matrixpower(x, power, complexOK = TRUE)
  matrixsqrt(x, complexOK = TRUE)
  matrixinvsqrt(x, complexOK = TRUE)
}
\arguments{
  \item{x}{
    A square matrix containing numeric or complex values.
  }
  \item{power}{
    A numeric value giving the power (exponent) to which \code{x} should
    be raised.
  }
  \item{complexOK}{
    Logical value indicating whether the result is allowed to be complex.
  }
}
\details{
  These functions raise the matrix \code{x} to the desired power:
  \code{matrixsqrt} takes the square root, \code{matrixinvsqrt} takes
  the inverse square root, and \code{matrixpower} takes the specified
  power of \code{x}.
  
  Up to numerical error, \code{matrixpower(x, 2)} should be equivalent
  to \code{x \%*\% x}, and \code{matrixpower(x, -1)} should be
  equivalent to \code{solve(x)}, the inverse of \code{x}. 

  The square root \code{y <- matrixsqrt(x)} should satisfy
  \code{y \%*\% y = x}. The inverse square root
  \code{z <- matrixinvsqrt(x)} should satisfy \code{z \%*\% z = solve(x)}.

  Computations are performed using the eigen decomposition
  (\code{\link{eigen}}).
}
\value{
  A matrix of the same size as \code{x} containing
  numeric or complex values.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[base]{eigen}}, \code{\link[base]{svd}}
}
\examples{
  x <- matrix(c(10,2,2,1), 2, 2)
  y <- matrixsqrt(x)
  y
  y \%*\% y 
  z <- matrixinvsqrt(x)
  z \%*\% y
  matrixpower(x, 0.1)
}
\keyword{algebra}
\keyword{array}
