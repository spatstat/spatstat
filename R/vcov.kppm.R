#
# vcov.kppm
#
#  vcov method for kppm objects
#
#   Original code: Abdollah Jalilian
#
#   $Revision: 1.7 $  $Date: 2015/02/20 11:54:41 $
#

vcov.kppm <- function(object, ...,
                      what=c("vcov", "corr", "fisher", "internals"),
                      fast = NULL, rmax = NULL, eps.rmax = 0.01,
                      verbose = TRUE)
{
  what <- match.arg(what)
  verifyclass(object, "kppm")
  if(is.null(object$improve)){#Normal composite likelihood (poisson) case
  # extract composite likelihood results
  po <- object$po
  # ensure it was fitted with quadscheme
  if(is.null(getglmfit(po))) {
    warning("Re-fitting model with forcefit=TRUE")
    po <- update(po, forcefit=TRUE)
  }
  # extract quadrature scheme information
  Q <- quad.ppm(po)
  U <- union.quad(Q)
  nU <- npoints(U)
  wt <- w.quad(Q)
  # compute fitted intensity values
  lambda <- fitted(po, type="lambda")
  # extract covariate values
  Z <- model.matrix(po)
  # extract pcf
  g <- pcfmodel(object)
  # resolve options for algorithm
  fast.given <- !is.null(fast)
  maxmat <- spatstat.options("maxmatrix")
  blocking <- FALSE
  if(!fast.given) {
    fast <- (nU^2 > maxmat)
  } else stopifnot(is.logical(fast))
  # decide whether matrix must be split into blocks
  if(!fast) {
    blocking <- (nU^2 > maxmat)
  } else {
    if(is.null(rmax)){
        diamwin <- diameter(as.owin(U))
        fnc <- get("fnc", envir = environment(improve.kppm))
        rmax <- if(fnc(diamwin, eps.rmax, g) >= 0) diamwin else
                  uniroot(fnc, lower = 0, upper = diamwin,
                          eps=eps.rmax, g=g)$root
    }
    ## approximate number of pairs
    nguess <- ceiling(4 * pi * (nU^2) * (rmax^2)/area(Window(U)))
    blocking <- (nguess > maxmat)
    if(blocking && fast.given)
      warning("Tapering algorithm could not be used: insufficient space.",
              " Reverting to full calculation.", call.=FALSE)
  }
  # evaluate integrand
  ff <- Z * lambda * wt
  ## compute quadratic form involving (pair correlation function minus 1)
  if(!blocking) {
    if(!fast) {
      gminus1 <- matrix(g(c(pairdist(U))) - 1, nU, nU)
    } else {
      cp <- crosspairs(U,U,rmax)
      gminus1 <- Matrix::sparseMatrix(i=cp$i, j=cp$j,
                                      x=g(cp$d) - 1, dims=c(nU, nU))
    }
    E <- t(ff) %*% gminus1 %*% ff
  } else {
    ## split calculation of (gminus1 %*% ff) into blocks
    nrowperblock <- max(1, floor(maxmat/nU))
    nblocks <- ceiling(nU/nrowperblock)
    g1ff <- NULL
    if(verbose)
      splat("Splitting large matrix calculation into", nblocks, "blocks")
    for(k in seq_len(nblocks)) {
      if(verbose) progressreport(k, nblocks)
      ii <- nrowperblock * (k-1) + seq_len(nrowperblock)
      ii <- ii[ii <= nU]
      gm1 <- matrix(g(c(crossdist(U[ii], U))) - 1, length(ii), nU)
      g1ff <- rbind(g1ff, gm1 %*% ff)
    }
    E <- t(ff) %*% g1ff
  }
  # asymptotic covariance matrix in the Poisson case
  J <- t(Z) %*% ff
  J.inv <- try(solve(J))
  # could be singular 
  if(inherits(J.inv, "try-error")) {
    if(what == "internals") {
      return(list(ff=ff, J=J, E=E, J.inv=NULL))
    } else {
      return(NULL)
    }
  }
  # asymptotic covariance matrix in the clustered case
  vc <- J.inv + J.inv %*% E %*% J.inv
  #
  } else {
    ## Case of quasi-likelihood (or other things from improve.kppm)
    run <- is.null(object$vcov) ||
      (!is.null(fast) && (fast != object$improve$fast.vcov))
    if(run){
      ## Calculate vcov if it hasn't already been so
      ## or if option fast differs from fast.vcov
      args <- object$improve
      internal <- what=="internals"
      if(!is.null(fast))
        args$fast.vcov <- fast
      object <- with(args,
                     improve.kppm(object, type = type,
                                  rmax = rmax, dimyx = dimyx,
                                  fast = fast, vcov = TRUE,
                                  fast.vcov = fast.vcov,
                                  maxIter = 0,
                                  save.internals = internal))
    }
    vc <- object$vcov
  }

  ## Convert from Matrix to ordinary matrix:
  vc <- as.matrix(vc)
  
  switch(what,
         vcov={ return(vc) },
         corr={
           sd <- sqrt(diag(vc))
           co <- vc/outer(sd, sd, "*")
           return(co)
         },
         fisher={
           fish <- try(solve(vc))
           if(inherits(fish, "try-error")) fish <- NULL 
           return(fish)
         },
         internals={
           return(list(ff=ff, J=J, E=E, J.inv=J.inv, vc=vc))
         })
  stop(paste("Unrecognised option: what=", what))
}
