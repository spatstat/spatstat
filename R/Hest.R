#
#    Hest.R
#
#  Contact distribution for a random set
#
#
Hest <- local({
  Hest <- function(X, r=NULL, breaks=NULL,
                   ...,
                   W,
                   correction=c("km", "rs", "han"),
                   conditional=TRUE) {
    rorbgiven <- !is.null(r) || !is.null(breaks)
    checkspacing <- !isFALSE(list(...)$checkspacing)
    testme       <- isTRUE(list(...)$testme)
    if(is.ppp(X) || is.psp(X)) {
      XX <- X
      W0 <- Window(X)
    } else if(is.owin(X)) {
      XX <- X
      W0 <- Frame(X)
    } else if(is.im(X)) {
      if(!is.logical(ZeroValue(X)))
        stop("When X is an image, its pixel values should be logical values")
      XX <- solutionset(X)
      W0 <- Window(X)
    } else stop("X should be an object of class ppp, psp, owin or im")
    ##
    if(given.W <- !missing(W) && !is.null(W)) {
      stopifnot(is.owin(W))
      if(!is.subset.owin(W, W0))
        stop("W is not a subset of the observation window of X")
    } else {
      W <- W0
    }
    ## handle corrections
    if(is.null(correction))
      correction <- c("rs", "km", "cs")
    correction <- pickoption("correction", correction,
                             c(none="none",
                               raw="none",
                               border="rs",
                               rs="rs",
                               KM="km",
                               km="km",
                               Kaplan="km",
                               han="han",
                               Hanisch="han",
                               best="km"),
                             multi=TRUE)
    corxtable <- c("km", "rs", "han", "none") 
    corx <- as.list(corxtable %in% correction)
    names(corx) <- corxtable
    ## compute distance map
    D <- distmap(XX, ...)
    pixeps <- with(D, min(xstep, ystep))
    if(!given.W && !is.im(X)) {
      B <- attr(D, "bdry")
    } else {
      B <- distmap(W, invert=TRUE, ...)
      har <- harmonise(D=D, B=B)
      D <- har$D[W, drop=FALSE]
      B <- har$B[W, drop=FALSE]
    }
    ## histogram breakpoints 
    dmax <- max(D)
    breaks <- handle.r.b.args(r, breaks, W, NULL, rmaxdefault=dmax)
    rval <- breaks$r
    if(testme || (rorbgiven && checkspacing))
      check.finespacing(rval, rname="r", eps=pixeps/4, W,
                        rmaxdefault=dmax,
                        context="in Hest(X,r)",
                        action="fatal")
    ##  extract distances and censoring distances
    dist <- as.vector(as.matrix(D))
    bdry <- as.vector(as.matrix(B))
    ok <- !is.na(dist) & !is.na(bdry)
    dist <- dist[ok]
    bdry <- bdry[ok]
    ## delete zero distances
    if(is.owin(X) || is.im(X)) {
      pos <- (dist > 0)
      areafraction <- 1 - mean(pos)
      dist <- dist[pos]
      bdry <- bdry[pos]
    }
    ## censoring indicators
    d <- (dist <= bdry)
    ##  observed distances
    o <- pmin.int(dist, bdry)
    ## calculate estimates
    Z <- censtimeCDFest(o, bdry, d, breaks,
                        KM=corx$km,
                        RS=corx$rs,
                        HAN=corx$han,
                        RAW=corx$none,
                        han.denom=if(corx$han) eroded.areas(W, rval) else NULL,
                        tt=dist)
    ## conditional on d > 0 ?
    if(is.owin(X) || is.im(X)) {
      if(conditional) {
        if(corx$km)   Z$km  <- condition(Z$km)
        if(corx$rs)   Z$rs  <- condition(Z$rs)
        if(corx$han)  Z$han <- condition(Z$han)
        if(corx$none) Z$raw <- condition(Z$raw)
      } else {
        if(corx$km)   Z$km  <- reconstitute(Z$km, areafraction) 
        if(corx$rs)   Z$rs  <- reconstitute(Z$rs, areafraction) 
        if(corx$han)  Z$han <- reconstitute(Z$han, areafraction) 
        if(corx$none) Z$raw <- reconstitute(Z$raw, areafraction) 
      }
    }
    ## relabel
    Z <- rebadge.fv(Z, substitute(H(r), NULL), "H")
    unitname(Z) <- unitname(X)
    attr(Z, "conserve") <- list(checkspacing=FALSE)
    return(Z)
  }

  condition <- function(x) { (x - x[1])/(1-x[1]) }
  reconstitute <- function(x, p) { p + (1-p) * x }

  Hest
})


