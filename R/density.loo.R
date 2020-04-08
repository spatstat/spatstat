#'
#'    density.loo.R
#'
#'  Compute leave-one-out density estimates at each data point
#'  on a linear network
#'
#'  Copyright (c) Greg McSwiggan and Adrian Baddeley 2017-2020

densitypointsLPP <- function(x, sigma, ..., 
                             weights=NULL, nsigma=1,
                             leaveoneout=TRUE, fast=TRUE,
                             fastmethod=c("onestep", "absorb"),
                             floored=TRUE,
                             dx=NULL, dt=NULL,
                             iterMax=1e6, verbose=FALSE, debug=FALSE) {
  stopifnot(is.lpp(x))
  #' compute density estimates at points
  fastmethod <- match.arg(fastmethod)
  if(identical(sigma, Inf)) {
    if(nsigma != 1) stop("nsigma should be equal to 1 when sigma is infinite")
    return(flatdensityatpointslpp(x, leaveoneout=leaveoneout,
                               weights=weights, disconnect=TRUE))
  }
  if(!leaveoneout || fast) {
    #' compute density estimates at points WITHOUT omitting points
    f <- densityfun(x, sigma, weights=weights, nsigma=nsigma, ..., 
                    dx=dx, dt=dt, iterMax=iterMax, verbose=verbose, debug=debug)
    tau <- attr(f, "sigma")
    y <- f(x)
    attr(y, "sigma") <- tau
    if(!leaveoneout)
      return(y)
    #' fast approximation to leave-one-out estimate
    #' evaluate approximation to heat kernel and subtract from 'y'
    if(nsigma == 1) y <- matrix(y, ncol=1)
    lenf <- lengths_psp(as.psp(domain(x)))
    coo <- coords(x)
    seg <- coo$seg
    len <- lenf[seg]
    pos <- len * coo$tp
    result <- ker0 <- matrix(, npoints(x), ncol(y))
    rr <- as.integer(row(ker0))
    cc <- as.integer(col(ker0))
    switch(fastmethod,
           absorb = {
             ker0[] <- hotrod(len[rr], pos[rr], pos[rr], tau[cc],
                              ends="absorbing")
           },
           onestep = {
             antipos <- len - pos
             L <- domain(x)
             vv <- vertexdegree(L)
             dleft <- vv[L$from[seg]]
             dright <- vv[L$to[seg]]
             Tau <- tau[cc]
             ker0[] <- pmax(0,
                            dnorm(0, 0, Tau) +
                            (2/dleft - 1) * dnorm(2 * pos[rr], 0, Tau) +
                            (2/dright - 1) * dnorm(2 * antipos[rr], 0, Tau))
           })
    if(floored)
      ker0[] <- pmax(ker0[], 1/volume(domain(x)))
    result[] <- pmax(0, y - ker0)
    if(ncol(result) == 1) result <- as.numeric(result)
    attr(result, "sigma") <- tau
    return(result)
  }
  #' get setup data
  g <- densityfun(x, sigma, weights=weights,
                  dx=dx, dt=dt, iterMax=iterMax,
                  verbose=verbose, debug=debug, 
                  exit="setup")
  #' extract internal data
  finenet     <- g$linnet_obj
  lixelmap    <- as.integer(g$lixelmap)
  lixelweight <- g$lixelweight
  Amatrix     <- g$Amatrix
  U0          <- g$U0
  deltax      <- g$deltax
  deltat      <- g$deltat
  #' 
  if(debug) {
    cat("finenet:\n")
    print(finenet)
    cat("lixelmap:\n")
    str(lixelmap)
    cat("lixelweight:\n")
    str(lixelweight)
  }
  #'
  niter <- ceiling(sigma^2/(deltat * 2))
  #' do the full iterative computation for each X[-i]
  U0 <- as.numeric(U0)
  v <- looHeatLPP(U0, Amatrix, npoints(x), niter, nsigma, lixelweight, lixelmap)
  result <- if(nsigma == 1) as.numeric(v) else t(v)
  attr(result, "sigma") <- sigma * sqrt(attr(v, "iter")/niter)
  return(result)
}


looHeatLPP <- function(U0, Amatrix, npts, niter, nsave, lixelweight, lixelmap,
                       verbose=TRUE) {
  pstate <- list()
  result <- matrix(nrow=nsave, ncol=npts)
  lixelmap <- as.integer(lixelmap)
  if(verbose) cat(paste("Processing", npts, "points ... "))
  if(nsave == niter) {
    #' save results of all iterations
    for(i in 1:npts) {
      u0 <- U0
      #' subtract weights corresponding to i-th data point
      ii <- 2 * i + c(-1, 0)
      ll <- lixelmap[ii]
      ww <- lixelweight[ii]
      u0[ll] <- u0[ll] - ww
      #' which node is closest to i-th data point
      kk <- ll[which.max(abs(ww))]
      #' run solver
      U <- u0
      for(iter in 1:niter) {
        U <- Amatrix %*% U
        result[iter, i] <- U[kk]
      }
      if(verbose) pstate <- progressreport(i, npts, state=pstate)
    }
    attr(result, "iter") <- 1:niter
  } else {
    #' save results of 'nsave' equally-spaced iterations
    blocksize <- ceiling(niter/nsave)
    for(i in 1:npts) {
      u0 <- U0
      #' subtract weights corresponding to i-th data point
      ii <- 2 * i + c(-1, 0)
      ll <- lixelmap[ii]
      ww <- lixelweight[ii]
      u0[ll] <- u0[ll] - ww
      #' which node is closest to i-th data point
      kk <- ll[which.max(abs(ww))]
      #' run solver
      U <- u0
      for(isave in 1:nsave) {
        nit <- min(blocksize, niter - (isave-1)*blocksize)
        for(j in 1:nit)
          U <- Amatrix %*% U
        result[isave, i] <- U[kk]
      }
      if(verbose) pstate <- progressreport(i, npts, state=pstate)
    }
    attr(result, "iter") <- pmin(niter, blocksize * (1:nsave))
  }
  return(result)
}
