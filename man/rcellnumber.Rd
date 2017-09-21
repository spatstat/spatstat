\name{rcellnumber}
\alias{rcellnumber}
\title{
  Generate Random Numbers of Points for Cell Process
}
\description{
  Generates random integers for the Baddeley-Silverman counterexample.
}
\usage{
 rcellnumber(n, N = 10, mu=1)
}
\arguments{
  \item{n}{
    Number of random integers to be generated.
  }
  \item{N}{
    Distributional parameter: the largest possible value
    (when \code{mu <= 1}).
    An integer greater than 1.
  }
  \item{mu}{
    Mean of the distribution (equals the variance).
    Any positive real number.
  }
}
\details{
  If \code{mu = 1} (the default),
  this function generates random integers which have mean and variance
  equal to 1, but which do not have a Poisson distribution. 
  The random integers take the values \eqn{0}, \eqn{1} and \eqn{N}
  with probabilities \eqn{1/N}, \eqn{(N-2)/(N-1)} and \eqn{1/(N(N-1))}
  respectively.
  See Baddeley and Silverman (1984).

  If \code{mu} is another positive number, the random integers will
  have mean and variance equal to \code{mu}. They are obtained by
  generating the
  one-dimensional counterpart of the cell process and counting the
  number of points in the interval from \code{0} to \code{mu}. The
  maximum possible value of each random integer is \code{N * ceiling(mu)}.
}
\value{
  An integer vector of length \code{n}.
}
\references{
  Baddeley, A.J. and Silverman, B.W. (1984)
  A cautionary example on the use of second-order methods for analyzing
  point patterns. \emph{Biometrics} \bold{40}, 1089-1094.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{rcell}}
}
\examples{
   rcellnumber(30, 3)
}
\keyword{datagen}
