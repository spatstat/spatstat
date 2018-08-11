#
#   rknn.R
#
#   Distribution of distance to k-th nearest point in d dimensions
#   (Poisson process of intensity lambda)
#
#   $Revision: 1.2 $  $Date: 2009/12/31 01:33:44 $
#

dknn <- function(x, k=1, d=2, lambda=1) {
  validposint(k, "dknn")
  validposint(d, "dknn")
  alpha.d <- (2 * pi^(d/2))/(d * gamma(d/2.))
  y <- dgamma(x^d, shape=k, rate=lambda * alpha.d)
  y <- y * d * x^(d-1)
  return(y)
}

pknn <- function(q, k=1, d=2, lambda=1) {
  validposint(k, "pknn")
  validposint(d, "pknn")
  alpha.d <- (2 * pi^(d/2))/(d * gamma(d/2.))
  p <- pgamma(q^d, shape=k, rate=lambda * alpha.d)
  return(p)
}

qknn <- function(p, k=1, d=2, lambda=1) {
  validposint(k, "qknn")
  validposint(d, "qknn")
  alpha.d <- (2 * pi^(d/2))/(d * gamma(d/2.))
  y <- qgamma(p, shape=k, rate=lambda * alpha.d)
  z <- y^(1/d)
  return(z)
}

rknn <- function(n, k=1, d=2, lambda=1) {
  validposint(k, "rknn")
  validposint(d, "rknn")
  alpha.d <- (2 * pi^(d/2))/(d * gamma(d/2.))  
  y <- rgamma(n, shape=k, rate=lambda * alpha.d)
  x <- y^(1/d)
  return(x)
}

  
