#
#   varblock.R
#
#   Variance estimation using block subdivision
#
#   $Revision: 1.8 $  $Date: 2010/11/21 07:05:12 $
#

varblock <- function(X, fun=Kest,
                     blocks=quadrats(X, nx=nx, ny=ny), ..., 
                     nx=3, ny=nx) {
  stopifnot(is.ppp(X))
  stopifnot(is.tess(blocks))
  if(is.character(fun)) fun <- get(fun, mode="function")
  stopifnot(is.function(fun))
  rvalues <- function(z) { with(z, .x) }
  canrestrict <- "domain" %in% names(formals(fun))
  if(!canrestrict) {
    # divide data into disjoint blocks
    Y <- split(X, blocks)
    n <- length(Y)
    if(n <= 1) stop("Need at least 2 blocks")
    # apply 'fun' to each block
    if(any(c("r", "breaks") %in% names(list(...)))) {
      # r vector specified
      fX <- fun(X, ...)
      z <- lapply(Y, fun, ...)
    } else {
      # need to ensure compatible fv objects
      z <- lapply(Y, fun, ...)
      rmaxes <- unlist(lapply(z, function(x){ max(rvalues(x)) }))
      smallest <- which.min(rmaxes)
      r <- rvalues(z[[smallest]])
      z <- lapply(Y, fun, ..., r=r)
      fX <- fun(X, ..., r=r)
    }
  } else {
    # use 'domain' argument of 'fun' to compute contributions from each tile
    B <- tiles(blocks)
    n <- length(B)
    dofun <- function(domain, fun, Xpp, ...) {
      fun(Xpp, ..., domain=domain)
    }
    if(any(c("r", "breaks") %in% names(list(...)))) {
      # r vector specified
      fX <- fun(X, ...)
      z <- lapply(B, dofun, ..., fun=fun, Xpp=X)
    } else {
      # need to ensure compatible fv objects
      z <- lapply(B, dofun, ..., fun=fun, Xpp=X)
      rmaxes <- unlist(lapply(z, function(x){ max(rvalues(x)) }))
      smallest <- which.min(rmaxes)
      r <- rvalues(z[[smallest]])
      z <- lapply(B, dofun, ..., fun=fun, Xpp=X, r=r)
      fX <- fun(X, ..., r=r)
    }
  }
  # find columns that are common to all estimates
  zzz <- reconcile.fv(append(list(fX), z))
  fX <- zzz[[1]]
  z <- zzz[-1]
  # get info
  ylab <- attr(z[[1]], "ylab")
  yexp <- attr(z[[1]], "yexp")
  fname <- attr(z[[1]], "fname")
  # sample mean
  m <- meanlistfv(z)
  # sample variance
  sqdev <- lapply(z, function(x,m){ eval.fv((x-m)^2) }, m=m)
  v <- meanlistfv(sqdev)
  v <- eval.fv(v * n/(n-1))
  v <- rebadge.fv(v, new.ylab=ylab, new.fname=fname, new.yexp=yexp)
  # sample standard deviation
  sd <- eval.fv(sqrt(v))
  sd <- rebadge.fv(sd, new.ylab=ylab, new.fname=fname, new.yexp=yexp)
  # upper and lower limits
  sem <- eval.fv(sd/sqrt(n))
  upper <- eval.fv(fX + 2 * sem)
  upper <- rebadge.fv(upper, new.ylab=ylab, new.fname=fname, new.yexp=yexp)
  lower <- eval.fv(fX - 2 * sem)
  lower <- rebadge.fv(lower, new.ylab=ylab, new.fname=fname, new.yexp=yexp)
  # tack together 
  m <- prefixfv(m, "mean", "sample mean of")
  v <- prefixfv(v, "var", "estimated variance of")
  sd <- prefixfv(sd, "sd", "estimated standard deviation of")
  upper <- prefixfv(upper, "hi", "upper CI limit for")
  lower <- prefixfv(lower, "lo", "lower CI limit for")
  out <- cbind(fX,m,v,sd,upper,lower)
  # restrict r domain
  ok <- apply(is.finite(as.matrix(as.data.frame(out))), 1, all)
  rmax <- max(rvalues(out)[ok])
  alim <- attr(out, "alim")
  attr(out, "alim") <- c(0, min(rmax, alim[2]))
  return(out)
}

meanlistfv <- function(z) {
  # compute sample mean of a list of fv objects
  if(!is.list(z) || !all(unlist(lapply(z, is.fv))))
    stop("z should be a list of fv objects")
  n <- length(z)
  if(n == 0) return(NULL)
  a <- a1 <- z[[1]]
  if(n > 1) {
    for(i in 2:n) { 
      b <- z[[i]]
      a <- eval.fv(a+b) 
    }
    a <- eval.fv(a/n)
    a <- rebadge.fv(a, new.ylab=attr(a1, "ylab"),
                       new.fname=attr(a1, "fname"))
  }
  return(a)
}

  
