##
## bw.diggle.R
##
## bandwidth selection rules bw.diggle and bw.scott (for density.ppp)
##
## $Revision: 1.1 $ $Date: 2014/05/20 09:30:49 $
##

bw.scott <- function(X) {
  stopifnot(is.ppp(X))
  n <- npoints(X)
  sdx <- sqrt(var(X$x))
  sdy <- sqrt(var(X$y))
  return(c(sdx, sdy) * n^(-1/6))
}

bw.diggle <- function(X, ..., correction="good", hmax=NULL) {
  stopifnot(is.ppp(X))
  # secret option for debugging
  mf <- function(..., method=c("C", "interpreted")) match.arg(method)
  method <- mf(...)
  #
  lambda <- npoints(X)/area(Window(X))
  r <- if(is.null(hmax)) NULL else seq(0, 4*hmax, length=512)
  K <- Kest(X, r=r, correction=correction)
  yname <- fvnames(K, ".y")
  K <- K[, c("r", yname)]
  ## check that K values can be passed to C code
  if(any(bad <- !is.finite(K[[yname]]))) {
    ## throw out bad values
    lastgood <- min(which(bad)) - 1
    if(lastgood < 2)
      stop("K function yields too many NA/NaN values")
    K <- K[1:lastgood, ]
  }
  rvals <- K$r
  # evaluation of M(r) requires K(2r)
  rmax2 <- max(rvals)/2
  if(!is.null(alim <- attr(K, "alim"))) rmax2 <- min(alim[2], rmax2)
  ok <- (rvals <= rmax2)
  switch(method,
         interpreted = {
           rvals <- rvals[ok]
           nr <- length(rvals)
           J <- numeric(nr)
           phi <- function(x,h) { 
             if(h <= 0) return(numeric(length(x)))
             y <- pmax.int(0, pmin.int(1, x/(2 * h)))
             4 * pi * h^2 * (acos(y) - y * sqrt(1 - y^2))
           }
           for(i in 1:nr) 
             J[i] <- stieltjes(phi, K, h=rvals[i])[[yname]]/(2 * pi)
         },
         C = {
           nr <- length(rvals)
           nrmax <- sum(ok)
           dK <- diff(K[[yname]])
           ndK <- length(dK)
           DUP <- spatstat.options("dupC")
           z <- .C("digberJ",
                   r=as.double(rvals),
                   dK=as.double(dK),
                   nr=as.integer(nr),
                   nrmax=as.integer(nrmax),
                   ndK=as.integer(ndK),
                   J=as.double(numeric(nrmax)),
                   DUP=DUP)
           J <- z$J
           rvals <- rvals[ok]
         })
  pir2 <- pi * rvals^2
  M <- (1/lambda - 2 * K[[yname]][ok])/pir2 + J/pir2^2
  # This calculation was for the uniform kernel on B(0,h)
  # Convert to standard deviation of (one-dimensional marginal) kernel
  sigma <- rvals/2
  result <- bw.optim(M, sigma,
                     creator="bw.diggle",
                     criterion="Berman-Diggle Cross-Validation",
                     J=J,
                     lambda=lambda)
  return(result)
}


