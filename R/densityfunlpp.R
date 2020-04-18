##    densityfunlpp.R
##    Method for 'densityfun' for lpp objects
##
##    Copyright (c) Greg McSwiggan and Adrian Baddeley 2017-2020
##
##    $Revision: 1.10 $  $Date: 2020/04/18 08:10:11 $

densityfun.lpp <- function(X, sigma, ...,
                           weights=NULL, nsigma=1, verbose=FALSE) {
  stopifnot(is.lpp(X))
  check.1.real(sigma)
  if(sigma == Inf) {
    if(nsigma != 1)
      stop("nsigma must be equal to 1 when sigma is infinite")
    return(flatdensityfunlpp(X, weights=weights, disconnect=TRUE))
  } else check.finite(sigma)
  if(!is.null(weights)) 
    check.nvector(weights, npoints(X))
  #' 
  L <- as.linnet(X)
  p <- resolve.heat.steps(sigma, L=L, ..., nsave=nsigma, verbose=verbose)

  #' internal argument
  exit <- resolve.1.default(list(exit="no"), list(...))
  exit <- match.arg(exit, c("no", "parameters", "setup"))
  if(exit == "parameters") return(p)
  setuponly <- (exit == "setup")

  #' call Greg's solver
  a <- FDMKERNEL(lppobj=X, weights=weights, 
                 dtx=p$dx, dtt=p$dt, M=p$niter, nsave=p$nsave,
                 stepnames=list(time="dt", space="dx"),
                 setuponly=setuponly, verbose=verbose)
  if(setuponly) return(resolve.defaults(a, p))
  #' 
  if(nsigma == 1) {
    #' return smoother with bandwidth sigma
    result <- a$kernel_fun
    attr(result, "sigma") <- sigma
  } else {
    #' return multiple smoothers with bandwidths sigma * (k-1)/nsigma
    #' for k = 1, ..., nsigma+1
    result <- a$progressfun
    attr(result, "sigma") <- a$tau
  }
  attr(result, "dx") <- a$deltax
  attr(result, "dt") <- a$deltat
  return(result)
}
