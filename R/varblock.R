#
#   varblock.R
#
#   Variance estimation using block subdivision
#
#   $Revision: 1.19 $  $Date: 2016/03/08 00:48:06 $
#

varblock <- local({

  getrvalues <- function(z) { with(z, .x) }

  stepsize <- function(z) { mean(diff(z)) } 
  
  dofun <- function(domain, fun, Xpp, ...) { fun(Xpp, ..., domain=domain) }

  varblock <- function(X, fun=Kest,
                       blocks=quadrats(X, nx=nx, ny=ny), ..., 
                       nx=3, ny=nx,
                       confidence=0.95) {
    stopifnot(is.ppp(X))
    stopifnot(is.tess(blocks))
    stopifnot(is.function(fun) || is.character(fun))
    if(is.character(fun)) 
      fun <- get(fun, mode="function")
    ## validate confidence level
    stopifnot(confidence > 0.5 && confidence < 1)
    alpha <- 1 - confidence
    probs <- c(alpha/2, 1-alpha/2)
    ## determine whether 'fun' has an argument called 'domain'
    canrestrict <- ("domain" %in% names(formals(fun))) ||
                   samefunction(fun, pcf) ||
                   samefunction(fun, Lest)
    ## check there's at least one point in each block
    Y <- split(X, blocks)
    nums <- sapply(Y, npoints)
    blockok <- (nums > 0)
    if(some.zeroes <- any(!blockok)) 
      warning("Some tiles contain no data: they are discarded")
    if(!canrestrict) {
      ## divide data into disjoint blocks
      if(some.zeroes)
        Y <- Y[blockok]
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
        rlist <- lapply(z, getrvalues)
        rmax <- min(sapply(rlist, max))
        rstep <- min(sapply(rlist, stepsize))
        r <- seq(0, rmax, by=rstep)
        z <- lapply(Y, fun, ..., r=r)
        fX <- fun(X, ..., r=r)
      }
    } else {
      ## use 'domain' argument of 'fun' to compute contributions from each tile
      B <- tiles(blocks)
      if(some.zeroes)
        B <- B[blockok]
      n <- length(B)
      if(any(c("r", "breaks") %in% names(list(...)))) {
        ## r vector specified
        fX <- fun(X, ...)
        z <- lapply(B, dofun, ..., fun=fun, Xpp=X)
      } else {
        ## need to ensure compatible fv objects
        z <- lapply(B, dofun, ..., fun=fun, Xpp=X)
        rlist <- lapply(z, getrvalues)
        rmax <- min(sapply(rlist, max))
        rstep <- min(sapply(rlist, stepsize))
        r <- seq(0, rmax, by=rstep)
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
    sqdev <- lapply(z, sqdev.fv, m=m)
    v <- meanlistfv(sqdev)
    v <- eval.fv(v * n/(n-1), dotonly=FALSE)
    ## sample standard deviation
    sd <- eval.fv(sqrt(v), dotonly=FALSE)
    ## upper and lower limits
    sem <- eval.fv(sd/sqrt(n), dotonly=FALSE)
    zcrit <- qnorm(probs)
    lower <- eval.fv(m + zcrit[1] * sem, dotonly=FALSE)
    upper <- eval.fv(m + zcrit[2] * sem, dotonly=FALSE)
    ## rebadge
    fva <- .Spatstat.FvAttrib
    fva <- fva[fva %in% names(attributes(fX))]
    attributes(m)[fva] <- attributes(v)[fva] <- attributes(sd)[fva] <- 
        attributes(upper)[fva] <- attributes(lower)[fva] <- attributes(fX)[fva]
    m <- prefixfv(m, "mean", "sample mean of", "bold(mean)~")
    v <- prefixfv(v, "var", "estimated variance of", "bold(var)~")
    sd <- prefixfv(sd, "sd", "estimated standard deviation of", "bold(sd)~")
    CItext <- paste(c("lower", "upper"),
                    paste0(100 * confidence, "%%"),
                    "CI limit for")
    lower <- prefixfv(lower, "lo", CItext[1], "bold(lo)~")
    upper <- prefixfv(upper, "hi", CItext[2], "bold(hi)~")
    ## tack together 
    out <- cbind(fX,m,v,sd,upper,lower)
    ## restrict r domain
    bad <- apply(!is.finite(as.matrix(as.data.frame(out))), 1, all)
    rmax <- max(getrvalues(out)[!bad])
    alim <- c(0, rmax)
    if(!canrestrict) alim <- intersect.ranges(attr(out, "alim"), alim)
    attr(out, "alim") <- alim
    ## sensible default plot formula
    ybase <- fvnames(fX, ".y")
    xname <- fvnames(fX, ".x")
    tname <- intersect("theo", fvnames(fX, "."))
    fvnames(out, ".y") <- yname <- paste0("mean", ybase)
    fvnames(out, ".s") <- snames <- paste0(c("lo", "hi"), ybase)
    fvnames(out, ".") <- c(yname, tname, snames)
    attr(out, "fmla") <- paste(". ~ ", xname)
    return(out)
  }
  
  sqdev.fv <- function(x,m){ eval.fv((x-m)^2, dotonly=FALSE) }
  
  varblock
})


meanlistfv <- local({

  getYmatrix <- function(x, yn=ynames) { as.matrix(as.data.frame(x)[,yn]) }

  meanlistfv <- function(z, ...) {
    ## compute sample mean of a list of fv objects
    if(!is.list(z) || !all(unlist(lapply(z, is.fv))))
      stop("z should be a list of fv objects")
    if(!do.call(compatible, unname(z)))
      stop("Objects are not compatible")
    result <- template <- z[[1]]
    ## extract each object's function values as a matrix
    ynames <- fvnames(template, "*")
    matlist <- unname(lapply(z, getYmatrix, yn=ynames))
    ## stack matrices into an array
    y <- do.call(abind, append(matlist, list(along=3)))
    ## take mean 
    ymean <- apply(y, 1:2, mean, ...)
    result[,ynames] <- ymean
    return(result)
  }

  meanlistfv
})


  
