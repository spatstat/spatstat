#
# vcov.kppm
#
#  vcov method for kppm objects
#
#   Original code: Abdollah Jalilian
#
#   $Revision: 1.4 $  $Date: 2014/08/27 09:51:02 $
#

vcov.kppm <- function(object, ...,
                      what=c("vcov", "corr", "fisher", "internals"),
                      fast = NULL, rmax = NULL, eps.rmax = 0.01)
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
  # resolve fast option
  if(is.null(fast)){
    fast <- (nU > sqrt(spatstat.options("maxmatrix")))
  }
  if(!is.logical(fast))
      stop("Argument fast must be NULL or logical.")
  # compute pair correlation function minus 1
  if(fast){
    if(!require(Matrix))
      stop(paste("Package Matrix must be installed in order for",
                 "the fast option to work."),
           call.=FALSE)
    if(is.null(rmax)){
        diamwin <- diameter(as.owin(U))
        fnc <- get("fnc", envir = environment(improve.kppm))
        rmax <- if(fnc(diamwin, eps.rmax, g) >= 0) diamwin else
                  uniroot(fnc, lower = 0, upper = diamwin,
                          eps=eps.rmax, g=g)$root
    }
    cp <- crosspairs(U,U,rmax)
    gminus1 <- Matrix::sparseMatrix(i=cp$i, j=cp$j,
                                    x=g(cp$d) - 1, dims=c(nU, nU))
  } else{
    gminus1 <- matrix(g(c(pairdist(U))) - 1, nU, nU)
  }
  # evaluate integral
  ff <- Z * lambda * wt
  J <- t(Z) %*% ff
  E <- t(ff) %*% gminus1 %*% ff
  # asymptotic covariance matrix in the Poisson case
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
  } else{#Case of quasi-likelihood (or other things from improve.kppm)
      run <- is.null(object$vcov) || (!is.null(fast) && (fast!=object$improve$fast.vcov))
      if(run){#Calculate vcov if it hasn't already been so or option fast differs from fast.vcov
          args <- object$improve
          internal <- what=="internals"
          if(!is.null(fast))
              args$fast.vcov <- fast
          object <- with(args, improve.kppm(object, type = type,
                                            rmax = rmax, dimyx = dimyx,
                                            fast = fast, vcov = TRUE, fast.vcov = fast.vcov,
                                            maxIter = 0, save.internals = internal))
      }
      vc <- object$vcov
  }
  
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
