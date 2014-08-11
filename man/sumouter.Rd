\name{sumouter}
\alias{sumouter}
\alias{quadform}
\title{Compute Quadratic Forms}
\description{
  Calculates certain quadratic forms of matrices.
}
\usage{
  sumouter(x, w=NULL)
  quadform(x, v)
}
\arguments{
  \item{x}{A matrix, whose rows are the vectors in the quadratic form.}
  \item{w}{Optional vector of weights}
  \item{v}{Matrix determining the quadratic form}
}
\value{
  A vector or matrix.
}
\details{
  The matrix \code{x} will be interpreted as a collection of row vectors.
  The command \code{sumouter} computes the sum of the outer
  products of these vectors, weighted by the entries of \code{w}:
  \deqn{
    M = \sum_i w_i x_i x_i^\top
  }{
    M = sum[i] (w[i] * outer(x[i,], x[i,]))
  }
  where the sum is over all rows of \code{x}
  (after removing any rows containing \code{NA} or other non-finite
  values).
  If \code{w} is missing, the weights will be taken as 1.
  The result is a \eqn{p \times p}{p * p} matrix where
  \code{p = ncol(x)}.
  
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

  # See what happens with NA's
  x[3,2] <- NA
  sumouter(x, w)
  quadform(x, v)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{array}
 
 
