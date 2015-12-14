\name{rknn}
\alias{dknn}
\alias{pknn}
\alias{qknn}
\alias{rknn}
\title{
  Theoretical Distribution of Nearest Neighbour Distance
}
\description{
  Density, distribution function, quantile function and random
  generation for the random distance to the \eqn{k}th nearest neighbour
  in a Poisson point process in \eqn{d} dimensions.
}
\usage{
dknn(x, k = 1, d = 2, lambda = 1)
pknn(q, k = 1, d = 2, lambda = 1)
qknn(p, k = 1, d = 2, lambda = 1)
rknn(n, k = 1, d = 2, lambda = 1)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations to be generated.}
  \item{k}{order of neighbour.}
  \item{d}{dimension of space.}
  \item{lambda}{intensity of Poisson point process.}
}
\details{
  In a Poisson point process in \eqn{d}-dimensional space, let
  the random variable \eqn{R} be
  the distance from a fixed point to the \eqn{k}-th nearest random point,
  or the distance from a random point to the
  \eqn{k}-th nearest other random point.

  Then \eqn{R^d} has a Gamma distribution with shape parameter \eqn{k}
  and rate \eqn{\lambda * \alpha}{lambda * alpha} where
  \eqn{\alpha}{alpha} is a constant (equal to the volume of the
  unit ball in \eqn{d}-dimensional space).
  See e.g. Cressie (1991, page 61).

  These functions support calculation and simulation for the
  distribution of \eqn{R}. 
}
\value{
  A numeric vector:
  \code{dknn} returns the probability density,
  \code{pknn} returns cumulative probabilities (distribution function),
  \code{qknn} returns quantiles,
  and \code{rknn} generates random deviates.
}
\references{
  Cressie, N.A.C. (1991)
  \emph{Statistics for spatial data}.
  John Wiley and Sons, 1991.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\examples{
  x <- seq(0, 5, length=20)
  densities <- dknn(x, k=3, d=2)
  cdfvalues <- pknn(x, k=3, d=2)
  randomvalues <- rknn(100, k=3, d=2)
  deciles <- qknn((1:9)/10, k=3, d=2)
}
\keyword{spatial}
\keyword{distribution}
