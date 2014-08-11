##
##  hermite.R
##
##  Gauss-Hermite quadrature
##
##  $Revision: 1.4 $  $Date: 2014/04/13 08:24:52 $
##

HermiteCoefs <- function(order) {
  ## compute coefficients of Hermite polynomial (unnormalised)
  x <- 1
  if(order > 0) 
    for(n in 1:order)
      x <- c(0, 2 * x) - c(((0:(n-1)) * x)[-1], 0, 0)
  return(x)
}

gauss.hermite <- function(f, mu=0, sd=1, ..., order=5) {
  stopifnot(is.function(f))
  stopifnot(length(mu) == 1)
  stopifnot(length(sd) == 1)
  ## Hermite polynomial coefficients (un-normalised)
  Hn <- HermiteCoefs(order)
  Hn1 <- HermiteCoefs(order-1)
  ## quadrature points
  x <- sort(Re(polyroot(Hn)))
  ## weights
  Hn1x <- matrix(Hn1, nrow=1) %*% t(outer(x, 0:(order-1), "^"))
  w <- 2^(order-1) * factorial(order) * sqrt(pi)/(order * Hn1x)^2
  ## adjust
  ww <- w/sqrt(pi)
  xx <- mu + sd * sqrt(2) * x
  ## compute
  ans <- 0
  for(i in seq_along(x))
    ans <- ans + ww[i] * f(xx[i], ...)
  return(ans)
}

dmixpois <- local({

  dpoisG <- function(x, ..., k, g) dpois(k, g(x))

  function(x, mu, sd, invlink=exp, GHorder=5) 
    gauss.hermite(dpoisG, mu=mu, sd=sd, g=invlink, k=x, order=GHorder)
})

pmixpois <- local({
  ppoisG <- function(x, ..., q, g, lot) ppois(q, g(x), lower.tail=lot)

  function(q, mu, sd, invlink=exp, lower.tail = TRUE, GHorder=5) 
    gauss.hermite(ppoisG, mu=mu, sd=sd, g=invlink, q=q, order=GHorder,
                 lot=lower.tail)
})
  
qmixpois <- function(p, mu, sd, invlink=exp, lower.tail = TRUE, GHorder=5) {
  ## guess upper limit
  ## Guess upper and lower limits
  pmin <- min(p, 1-p)/2
  lam.hi <- invlink(qnorm(pmin, mean=max(mu), sd=max(sd), lower.tail=FALSE))
  lam.lo <- invlink(qnorm(pmin, mean=min(mu), sd=max(sd), lower.tail=TRUE))
  kmin <- qpois(pmin, lam.lo, lower.tail=TRUE)
  kmax <- qpois(pmin, lam.hi, lower.tail=FALSE)
  kk <- kmin:kmax
  pp <- pmixpois(kk, mu, sd, invlink, lower.tail=TRUE, GHorder)
  ans <- if(lower.tail) kk[findInterval(p, pp, all.inside=TRUE)] else
         rev(kk)[findInterval(1-p, rev(1-pp), all.inside=TRUE)]
  return(ans)
}
  
rmixpois <- function(n, mu, sd, invlink=exp) {
  lam <- invlink(rnorm(n, mean=mu, sd=sd))
  y <- rpois(n, lam)
  return(y)
}

