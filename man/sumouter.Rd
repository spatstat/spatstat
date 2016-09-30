\name{sumouter}
\alias{sumouter}
\alias{quadform}
\alias{bilinearform}
\title{Compute Quadratic Forms}
\description{
  Calculates certain quadratic forms of matrices.
}
\usage{
  sumouter(x, w=NULL, y=x)
  quadform(x, v)
  bilinearform(x, v, y)
}
\arguments{
  \item{x,y}{A matrix, whose rows are the vectors in the quadratic form.}
  \item{w}{Optional vector of weights}
  \item{v}{Matrix determining the quadratic form}
}
\value{
  A vector or matrix.
}
\details{
  The matrices \code{x} and \code{y} will be interpreted as
  collections of row vectors. They must have the same number of rows.

  The command \code{sumouter} computes the sum of the outer
  products of corresponding row vectors, weighted by the entries of \code{w}:
  \deqn{
    M = \sum_i w_i x_i y_i^\top
  }{
    M = sum[i] (w[i] * outer(x[i,], y[i,]))
  }
  where the sum is over all rows of \code{x} 
  (after removing any rows containing \code{NA} or other non-finite
  values).
  If \code{w} is missing, the weights will be taken as 1.
  The result is a \eqn{p \times q}{p * q} matrix where
  \code{p = ncol(x)} and \code{q = ncol(y)}.
  
  The command \code{quadform} evaluates the quadratic form, defined by
  the matrix \code{v}, for each of the row vectors of \code{x}:
 \deqn{ 
   y_i = x_i V x_i^\top
  }{
    y[i] = x[i,] \%*\% v \%*\% t(x[i,])
  }
  The result \code{y} is a numeric vector of length \code{n} where
  \code{n = nrow(x)}. If \code{x[i,]} contains \code{NA} or
  other non-finite values, then \code{y[i] = NA}.

  The command \code{bilinearform} evaluates the more general bilinear
  form defined by the matrix \code{v}. Here \code{x} and \code{y} must
  be matrices of the same dimensions. For each row vector of
  \code{x} and corresponding row vector of \code{y}, the bilinear form is
 \deqn{ 
   z_i = x_i V y_i^\top
  }{
    z[i] = x[i,] \%*\% v \%*\% t(y[i,])
  }
  The result \code{z} is a numeric vector of length \code{n} where
  \code{n = nrow(x)}. If \code{x[i,]} or \code{y[i,]} contains \code{NA} or
  other non-finite values, then \code{z[i] = NA}.
}
\examples{
  x <- matrix(1:12, 4, 3)
  dimnames(x) <- list(c("Wilma", "Fred", "Barney", "Betty"), letters[1:3])
  x

  sumouter(x)

  w <- 4:1
  sumouter(x, w)
  v <- matrix(1, 3, 3)
  quadform(x, v)

  # should be the same as quadform(x, v)
  bilinearform(x, v, x)

  # See what happens with NA's
  x[3,2] <- NA
  sumouter(x, w)
  quadform(x, v)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{array}
 
 
