##
## bw.diggle.R
##
## bandwidth selection rules bw.diggle and bw.scott (for density.ppp)
##
## $Revision: 1.4 $ $Date: 2015/04/02 02:01:01 $
##

bw.scott <- function(X) {
  stopifnot(is.ppp(X))
  n <- npoints(X)
  sdx <- sqrt(var(X$x))
  sdy <- sqrt(var(X$y))
  return(c(sdx, sdy) * n^(-1/6))
}

bw.diggle <- local({

  #' integrand 
  phi <- function(x,h) { 
    if(h <= 0) return(numeric(length(x)))
    y <- pmax.int(0, pmin.int(1, x/(2 * h)))
    4 * pi * h^2 * (acos(y) - y * sqrt(1 - y^2))
  }
  
  #' secret option for debugging
  mf <- function(..., method=c("C", "interpreted")) { match.arg(method) }

  
  bw.diggle <- function(X, ..., correction="good", hmax=NULL, nr=512) {
    stopifnot(is.ppp(X))
    method <- mf(...)
    W <- Window(X)
    lambda <- npoints(X)/area(W)
    rmax <- if(!is.null(hmax)) (4 * hmax) else rmax.rule("K", W, lambda)
    r <- seq(0, rmax, length=nr)
    K <- Kest(X, r=r, correction=correction)
    yname <- fvnames(K, ".y")
    K <- K[, c("r", yname)]
    ## check that K values can be passed to C code
    if(any(bad <- !is.finite(K[[yname]]))) {
      ## throw out bad values
      lastgood <- min(which(bad)) - 1L
      if(lastgood < 2L)
        stop("K function yields too many NA/NaN values")
      K <- K[1:lastgood, ]
    }
    rvals <- K$r
    ## evaluation of M(r) requires K(2r)
    rmax2 <- max(rvals)/2
    if(!is.null(alim <- attr(K, "alim"))) rmax2 <- min(alim[2L], rmax2)
    ok <- (rvals <= rmax2)
    switch(method,
           interpreted = {
             rvals <- rvals[ok]
             nr <- length(rvals)
             J <- numeric(nr)
             for(i in 1:nr) 
               J[i] <- stieltjes(phi, K, h=rvals[i])[[yname]]/(2 * pi)
           },
           C = {
             nr <- length(rvals)
             nrmax <- sum(ok)
             dK <- diff(K[[yname]])
             ndK <- length(dK)
             z <- .C("digberJ",
                     r=as.double(rvals),
                     dK=as.double(dK),
                     nr=as.integer(nr),
                     nrmax=as.integer(nrmax),
                     ndK=as.integer(ndK),
                     J=as.double(numeric(nrmax)),
                     PACKAGE = "spatstat")
             J <- z$J
             rvals <- rvals[ok]
           })
    pir2 <- pi * rvals^2
    M <- (1/lambda - 2 * K[[yname]][ok])/pir2 + J/pir2^2
    ## This calculation was for the uniform kernel on B(0,h)
    ## Convert to standard deviation of (one-dimensional marginal) kernel
    sigma <- rvals/2
    result <- bw.optim(M, sigma,
                       creator="bw.diggle",
                       criterion="Berman-Diggle Cross-Validation",
                       J=J,
                       lambda=lambda,
                       unitname=unitname(X))
    return(result)
  }

  bw.diggle
})



