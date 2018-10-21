#
#  envelopelpp.R
#
#  $Revision: 1.24 $   $Date: 2018/10/21 09:58:37 $
#
#  Envelopes for 'lpp' objects
#
#

envelope.lpp <-
  function(Y, fun=linearK, nsim=99, nrank=1, ...,
           funargs=list(), funYargs=funargs,
           simulate=NULL, fix.n=FALSE, fix.marks=FALSE,
           verbose=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL,
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE,
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2,
           Yname=NULL, do.pwrong=FALSE, envir.simul=NULL) {
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- linearK

  if("clipdata" %in% names(list(...)))
    stop(paste("The argument", sQuote("clipdata"),
               "is not available for envelope.lpp"))
  
  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())
  
  if(!is.null(simulate)) {
    # ...................................................
    # Simulations are determined by 'simulate' argument
    # Processing is deferred to envelopeEngine
    simrecipe <- simulate
    # Data pattern is argument Y
    X <- Y
  } else if(!fix.n && !fix.marks) {
    # ...................................................
    # Realisations of complete spatial randomness
    # Data pattern X is argument Y
    # Data pattern determines intensity of Poisson process
    X <- Y
    nY <- npoints(Y)
    Yintens <- intensity(unmark(Y))
    Ymarx <- marks(Y)
    NETWORK <- Y$domain
    dont.complain.about(nY, Yintens, NETWORK)
    ## expression that will be evaluated
    simexpr <- if(is.null(Ymarx)) {
      #' unmarked point pattern
      expression(rpoislpp(Yintens, NETWORK))
    } else if(is.null(dim(Ymarx))) {
      #' single column of marks
      expression({
        A <- rpoislpp(Yintens, NETWORK);
        j <- sample(nY, npoints(A), replace=TRUE);
        A %mark% Ymarx[j]
      })
    } else {
      #' multiple columns of marks
      expression({
        A <- rpoislpp(Yintens, NETWORK);
        j <- sample(nY, npoints(A), replace=TRUE);
        A %mark% Ymarx[j, , drop=FALSE]
      })
    }
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type = "csr",
                             expr = simexpr,
                             envir = envir.here,
                             csr   = TRUE)
  } else if(!fix.marks) {
    # Fixed number of points, but random locations and marks
    # Data pattern X is argument Y
    X <- Y
    nY <- npoints(Y)
    Ymarx <- marks(Y)
    NETWORK <- Y$domain
    dont.complain.about(nY, NETWORK)
    # expression that will be evaluated
    simexpr <- if(is.null(Ymarx)) {
      ## unmarked
      expression(runiflpp(nY, NETWORK))
    } else if(is.null(dim(Ymarx))) {
      ## single column of marks
      expression({
        A <- runiflpp(nY, NETWORK);
        j <- sample(nY, npoints(A), replace=TRUE);
        A %mark% Ymarx[j]
      })
    } else {
      ## multiple columns of marks
      expression({
        A <- runiflpp(nY, NETWORK);
        j <- sample(nY, npoints(A), replace=TRUE);
        A %mark% Ymarx[j, ,drop=FALSE]
      })
    }
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type = "csr",
                             expr = simexpr,
                             envir = envir.here,
                             csr   = TRUE)
  } else {
    # ...................................................
    # Randomised locations only; 
    # fixed number of points and fixed marks
    # Data pattern X is argument Y
    X <- Y
    nY <- npoints(Y)
    Ymarx <- marks(Y)
    NETWORK <- Y$domain
    # expression that will be evaluated
    simexpr <- expression(runiflpp(nY, NETWORK) %mark% Ymarx)
    dont.complain.about(nY, Ymarx, NETWORK)
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type = "csr",
                             expr = simexpr,
                             envir = envir.here,
                             csr   = TRUE)
  }
  envelopeEngine(X=X, fun=fun, simul=simrecipe,
                 nsim=nsim, nrank=nrank, ...,
                 funargs=funargs, funYargs=funYargs,
                 verbose=verbose, clipdata=FALSE,
                 transform=transform,
                 global=global, ginterval=ginterval, use.theory=use.theory,
                 alternative=alternative, scale=scale, clamp=clamp, 
                 savefuns=savefuns, savepatterns=savepatterns, nsim2=nsim2,
                 VARIANCE=VARIANCE, nSD=nSD,
                 Yname=Yname, cl=cl,
                 envir.user=envir.user, do.pwrong=do.pwrong)
}

envelope.lppm <-
  function(Y, fun=linearK, nsim=99, nrank=1, ..., 
           funargs=list(), funYargs=funargs,
           simulate=NULL, fix.n=FALSE, fix.marks=FALSE, verbose=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL,
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE, 
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2,
           Yname=NULL, do.pwrong=FALSE, envir.simul=NULL) {
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- linearK

  if("clipdata" %in% names(list(...)))
    stop(paste("The argument", sQuote("clipdata"),
               "is not available for envelope.pp3"))

  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())
  
  if(!is.null(simulate)) {
    # ...................................................
    # Simulations are determined by 'simulate' argument
    # Processing is deferred to envelopeEngine
    simrecipe <- simulate
    # Data pattern is argument Y
    X <- Y
  } else {
    ## ...................................................
    ## Simulation of the fitted model Y
    if(!is.poisson(Y))
      stop("Simulation of non-Poisson models is not yet implemented")
    MODEL <- Y
    X <- Y$X
    NETWORK <- domain(X)
    lambdaFit <- predict(MODEL)
    Xmarx <- marks(X)
    nX <- if(!is.marked(X)) npoints(X) else table(marks(X))
    dont.complain.about(NETWORK, Xmarx, nX)
    #' 
    if(!fix.n && !fix.marks) {
      #' Unconstrained simulations
      LMAX <-
        if(is.im(lambdaFit)) max(lambdaFit) else sapply(lambdaFit, max)
      dont.complain.about(LMAX)
      simexpr <- expression(rpoislpp(lambdaFit, NETWORK, lmax=LMAX))
    } else if(!fix.marks && is.marked(X)) {
      #' Fixed total number of points
      EN <- sapply(lambdaFit, integral)
      PROB <- EN/sum(EN)
      dont.complain.about(PROB)
      simexpr <- expression(
        rlpp(as.integer(rmultinom(1L, nX, PROB)), lambdaFit)
      )
    } else {
      #' Fixed number of points of each type 
      simexpr <- expression(rlpp(nX, lambdaFit))
    }
    #' evaluate in THIS environment
    simrecipe <- simulrecipe(type = "lppm",
                             expr = simexpr,
                             envir = envir.here,
                             csr   = FALSE)
  }
  envelopeEngine(X=X, fun=fun, simul=simrecipe,
                 nsim=nsim, nrank=nrank, ...,
                 funargs=funargs, funYargs=funYargs,
                 verbose=verbose, clipdata=FALSE,
                 transform=transform,
                 global=global, ginterval=ginterval, use.theory=use.theory,
                 alternative=alternative, scale=scale, clamp=clamp,
                 savefuns=savefuns, savepatterns=savepatterns, nsim2=nsim2,
                 VARIANCE=VARIANCE, nSD=nSD,
                 Yname=Yname, cl=cl,
                 envir.user=envir.user, do.pwrong=do.pwrong)
}
