#'
#' rhomono.R
#'
#' monotonic regression of (intensity of) point pattern X on covariate Z
#'
#' Copyright (c) Adrian Baddeley 2018
#' GNU Public Licence >= 2.0
#'
#'  $Revision: 1.3 $ $Date: 2019/08/26 05:03:33 $

rhohatMonotone <- function(X, Z, increasing=FALSE) {
  stopifnot(is.ppp(X) || is.lpp(X))
  nullmodel <- if(is.ppp(X)) ppm(X) else lppm(X) 
  stuff <- evalCovar(nullmodel, Z)$values
  #' observed values at data points
  x <- sort(stuff$ZX)
  #' unnormalised cdf of covariate
  w <- stuff$weights
  z <- stuff$Zvalues
  if(length(w) == 1) w <- rep(w, length(z))
  g <- ewcdf(z, w, normalise=FALSE)
  #'
  areas <- g(x)
  if(increasing) areas <- sum(w) - rev(areas)
  ## maximum upper sets algorithm
  y <- numeric(0)
  b <- diff(c(0, areas))
  while(length(b) > 0) {
    u <- seq_along(b)/cumsum(b)
    if(any(bad <- !is.finite(u))) # divide by zero etc
      u[bad] <- max(u[!bad], 0)
    k <- which.max(u)
    y <- c(y, rep(u[k], k))
    b <- b[-(1:k)]
  }
  y <- c(y, 0)
  if(increasing) y <- rev(y)
  lambda <- stepfun(x = x, y=y, right=TRUE, f=1)
  attr(lambda, "call") <- sys.call()
  return(lambda)
}
