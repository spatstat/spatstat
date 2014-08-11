#
#   varblock.R
#
#   Variance estimation using block subdivision
#
#   $Revision: 1.11 $  $Date: 2014/02/15 04:08:01 $
#

varblock <- local({

  rvalues <- function(z) { with(z, .x) }

  dofun <- function(domain, fun, Xpp, ...) { fun(Xpp, ..., domain=domain) }

  varblock <- function(X, fun=Kest,
                       blocks=quadrats(X, nx=nx, ny=ny), ..., 
                       nx=3, ny=nx) {
    stopifnot(is.ppp(X))
    stopifnot(is.tess(blocks))
    stopifnot(is.function(fun) || is.character(fun))
    if(is.character(fun)) 
      fun <- get(fun, mode="function")
    ## need method rather than generic, to detect 'domain' argument
    if(samefunction(fun, pcf))
      fun <- pcf.ppp
    canrestrict <- "domain" %in% names(formals(fun))
    if(!canrestrict) {
      ## divide data into disjoint blocks
      Y <- split(X, blocks)
      n <- length(Y)
      if(n <= 1) stop("Need at least 2 blocks")
      ## apply 'fun' to each block
      if(any(c("r", "breaks") %in% names(list(...)))) {
        ## r vector specified
        fX <- fun(X, ...)
        z <- lapply(Y, fun, ...)
      } else {
        ## need to ensure compatible fv objects
        z <- lapply(Y, fun, ...)
        rmaxes <- unlist(lapply(z, function(x){ max(rvalues(x)) }))
        smallest <- which.min(rmaxes)
        r <- rvalues(z[[smallest]])
        z <- lapply(Y, fun, ..., r=r)
        fX <- fun(X, ..., r=r)
      }
    } else {
      ## use 'domain' argument of 'fun' to compute contributions from each tile
      B <- tiles(blocks)
      n <- length(B)
      if(any(c("r", "breaks") %in% names(list(...)))) {
        ## r vector specified
        fX <- fun(X, ...)
        z <- lapply(B, dofun, ..., fun=fun, Xpp=X)
      } else {
        ## need to ensure compatible fv objects
        z <- lapply(B, dofun, ..., fun=fun, Xpp=X)
        rmaxes <- unlist(lapply(z, function(x){ max(rvalues(x)) }))
        smallest <- which.min(rmaxes)
        r <- rvalues(z[[smallest]])
        z <- lapply(B, dofun, ..., fun=fun, Xpp=X, r=r)
        fX <- fun(X, ..., r=r)
      }
    }
    ## find columns that are common to all estimates
    zzz <- reconcile.fv(append(list(fX), z))
    fX <- zzz[[1]]
    z <- zzz[-1]
    ## sample mean
    m <- meanlistfv(z)
    ## sample variance
    sqdev <- lapply(z, function(x,m){ eval.fv((x-m)^2) }, m=m)
    v <- meanlistfv(sqdev)
    v <- eval.fv(v * n/(n-1))
    ## sample standard deviation
    sd <- eval.fv(sqrt(v))
    ## upper and lower limits
    sem <- eval.fv(sd/sqrt(n))
    upper <- eval.fv(fX + 2 * sem)
    lower <- eval.fv(fX - 2 * sem)
    ## rebadge
    attributes(m) <- attributes(v) <- attributes(sd) <- attributes(fX)
    attributes(upper) <- attributes(lower) <- attributes(fX)
    m <- prefixfv(m, "mean", "sample mean of", "bold(mean)~")
    v <- prefixfv(v, "var", "estimated variance of", "bold(var)~")
    sd <- prefixfv(sd, "sd", "estimated standard deviation of", "bold(sd)~")
    upper <- prefixfv(upper, "hi", "upper CI limit for", "bold(hi)~")
    lower <- prefixfv(lower, "lo", "lower CI limit for", "bold(lo)~")
    ## tack together 
    out <- cbind(fX,m,v,sd,upper,lower)
    ## restrict r domain
    ok <- apply(is.finite(as.matrix(as.data.frame(out))), 1, all)
    rmax <- max(rvalues(out)[ok])
    alim <- attr(out, "alim")
    attr(out, "alim") <- c(0, min(rmax, alim[2]))
    return(out)
  }

  varblock
})


meanlistfv <- function(z) {
  # compute sample mean of a list of fv objects
  if(!is.list(z) || !all(unlist(lapply(z, is.fv))))
    stop("z should be a list of fv objects")
  if(!do.call("compatible", unname(z)))
    stop("Objects are not compatible")
  result <- template <- z[[1]]
  ynames <- fvnames(template, "*")
  f <- function(x, ynames=ynames) { as.matrix(as.data.frame(x))[,ynames] }
  y <- array(unlist(lapply(z, f)),
             dim=c(nrow(result), ncol(result)-1, length(z)))
  ymean <- apply(y, 1:2, mean)
  result[,ynames] <- ymean
  return(result)
}

  
