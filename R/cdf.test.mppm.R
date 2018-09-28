#
# cdf.test.mppm.R
#
# $Revision: 1.17 $  $Date: 2018/09/28 05:13:38 $
#
cdf.test.mppm <- local({

  allpixelvalues <- function(z) { as.vector(as.matrix(z)) }

  xcoord <- function(x, y) { x }
  ycoord <- function(x, y) { y }
  
  cdf.test.mppm <- function(model, covariate,
                            test=c("ks", "cvm", "ad"), ...,
                            nsim=19, verbose=TRUE,
                            interpolate=FALSE, fast=TRUE, jitter=TRUE) {
    modelname <- short.deparse(substitute(model))
    covname <- short.deparse(substitute(covariate))
    test <- match.arg(test)
    result <- PoissonTest(model, covariate, test=test, ...,
                          verbose=FALSE,
                          interpolate=interpolate, fast=fast, jitter=jitter,
                          modelname=modelname, covname=covname,
                          gibbsok=TRUE)
    if(is.poisson(model))
      return(result)
    result$poisson.p.value <- pobs <- result$p.value
    result$poisson.statistic <- tobs <- result$statistic
    ## Simulate ...
    Sims <- simulate(model, nsim=nsim, ..., verbose=verbose)
    if(verbose) 
      cat("Processing ...")
    state <- list()
    Yname <- model$Info$Yname
    Data <- eval(getCall(model)$data,
                 envir=environment(terms(model)))
    sim.pvals <- sim.stats <- numeric(nsim)    
    for(isim in 1:nsim) {
      Data[,Yname] <- Sims[,isim,drop=FALSE]
      modeli <- update(model, data=Data)
      Ai <- PoissonTest(modeli, covariate, test=test, ...,
                        verbose=FALSE,
                        interpolate=interpolate, fast=fast, jitter=jitter,
                        modelname=modelname, covname=covname,
                        gibbsok=TRUE)
      sim.pvals[isim] <- Ai$p.value
      sim.stats[isim] <- Ai$statistic
      if(verbose) state <- progressreport(isim, nsim, state=state)
    }
    ### COMPUTE p-value and pack up
    result$sim.pvals <- sim.pvals
    result$sim.stats <- sim.stats
    ## Monte Carlo p-value
    ## For tied p-values, first compare values of test statistics
    ## (because p = 0 may occur due to rounding)
    ## otherwise resolve ties by randomisation
    nless <- sum(sim.pvals < pobs)
    nplus <- sum(sim.pvals == pobs & sim.stats > tobs)
    nties <- sum(sim.pvals == pobs & sim.stats == tobs) 
    result$p.value <- (nless + nplus + sample(0:nties, 1L))/(nsim+1L)
    ##
    result$method <- c("Monte Carlo test of fitted Gibbs model",
                       paste("based on", nsim, "repetitions of"),
                       sub("Spatial", "spatial", result$method))
    return(result)
  }

  PoissonTest <- function(model, covariate,
                          test=c("ks", "cvm", "ad"), ..., verbose=TRUE,
                          interpolate=FALSE, fast=TRUE, jitter=TRUE,
                          gibbsok=FALSE,
                          modelname, covname) {
    if(missing(modelname)) modelname <- short.deparse(substitute(model))
    if(missing(covname)) covname <- short.deparse(substitute(covariate))
    test <- match.arg(test)
    stopifnot(is.mppm(model))
    if(!gibbsok && !is.poisson.mppm(model))
      stop("Only implemented for Poisson models")
    ## extract things from model
    data  <- model$data
    npat  <- model$npat
    Y     <- data.mppm(model)
    if(fast) {
      ## extract original quadrature schemes and convert to point patterns
      QQ  <- quad.mppm(model)
      PP  <- lapply(QQ, union.quad)
      Zweights <- lapply(QQ, w.quad)
    } else
      Zweights <- list()
    ## `evaluate' covariate
    if(verbose)
      cat("Extracting covariate...")
    if(identical(covariate, "x")) covariate <- xcoord
    if(identical(covariate, "y")) covariate <- ycoord
    if(is.character(covariate)) {
      ## extract covariate with this name from data used to fit model
      if(!(covariate %in% names(data)))
        stop(paste("Model does not contain a covariate called",
                   dQuote(covariate)))
      covname <- covariate
      covariate <- data[, covname, drop=TRUE]
    } else if(inherits(covariate, c("listof", "anylist"))) {
      if(length(covariate) != npat)
        stop(paste("Length of list of covariate values does not match",
                   "number of point patterns in data of original model"))
    } else if(is.hyperframe(covariate)) {
      ## extract first column
      covariate <- covariate[,1L, drop=TRUE]
      if(length(covariate) != npat)
        stop(paste("Number of rows of covariate hyperframe does not match",
                   "number of point patterns in data of original model"))
    } else if(is.function(covariate) || is.im(covariate)) {
      ## replicate to make a list
      covariate <- as.anylist(rep(list(covariate), npat))
    } else     
      stop(paste("Format of argument", sQuote("covariates"), "not understood"))
    if(verbose) {
      cat("done.\nComputing statistics for each pattern...")
      pstate <- list()
    }
    ## compile information for test from each row
    Zvalues <- ZX <- Win <- list()
    for(i in 1:npat) {
      if(verbose) pstate <- progressreport(i, npat, state=pstate)
      XI <- Y[[i]]
      if(fast)
        PI <- PP[[i]]
      else
        WI <- XI$window
      covariateI <- covariate[[i]]
      if(is.im(covariateI)) {
        type <- "im"
        ## evaluate at data points
        ZXI <-
          if(interpolate) interp.im(covariateI, XI$x, XI$y)
          else covariateI[XI]
        if(fast) {
          ## covariate values for quadrature points
          ZI <- covariateI[PI]
        } else {
          ## covariate image inside window
          ZI <- covariateI[WI, drop=FALSE]
          ## corresponding mask
          WI <- as.owin(ZI)
          ## pixel areas 
          Zweights[[i]] <- rep(WI$xstep * WI$ystep, prod(WI$dim))
        }
      } else if(is.function(covariateI)) {
        type <- "function"
        ## evaluate exactly at data points
        ZXI <- covariateI(XI$x, XI$y)
        if(fast) {
          ## covariate values for quadrature points
          ZI <- covariateI(PI$x, PI$y)
        } else {
          ## window
          WI <- as.mask(WI)
          ## covariate image inside window
          ZI <- as.im(covariateI, W=WI)
          ## pixel areas 
          Zweights[[i]] <- rep(WI$xstep * WI$ystep, prod(WI$dim))
        }
      } else
        stop("covariate should be an image or a function(x,y)")
      ZX[[i]] <- ZXI
      if(fast)
        Zvalues[[i]] <- ZI      
      else {
        Win[[i]] <- WI
        ## values of covariate in window
        Zvalues[[i]] <- allpixelvalues(ZI)
      }
    }

    if(verbose)
      cat("done.\nComputing predicted intensity...")

    ## compute predicted intensities
    trend <-
      if(fast)
        fitted(model, type="trend")
      else
        predict(model, type="trend", locations=Win, verbose=verbose)$trend
  
    if(verbose)
      cat("done.\nExtracting...")
    ## extract relevant values
    lambda <- if(fast) trend else lapply(trend, allpixelvalues)
    if(verbose)
      cat("done.\nPerforming test...")
  
    ## flatten to vectors
    lambda <- unlist(lambda)
    Zweights <- unlist(Zweights)
    Zvalues <- unlist(Zvalues)
    ZX      <- unlist(ZX)
    if(length(lambda) != length(Zvalues))
      stop("Internal error: mismatch between predicted values and Z values")
    if(length(Zvalues) != length(Zweights))
      stop("Internal error: mismatch between Z values and Z weights")
    lambda <- lambda * Zweights
  
    ## form weighted cdf of Z values in window
    FZ <- ewcdf(Zvalues, lambda/sum(lambda))
    ## Ensure support of cdf includes the range of the data
    xxx <- knots(FZ)
    yyy <- FZ(xxx)
    if(min(xxx) > min(ZX)) {
      xxx <- c(min(ZX), xxx)
      yyy <- c(0, yyy)
    }
    if(max(xxx) < max(ZX)) {
      xxx <- c(xxx, max(ZX))
      yyy <- c(yyy, 1)
    }
    ## make piecewise linear approximation of cdf
    FZ <- approxfun(xxx, yyy, rule=2)
    ## evaluate at data points
    if(!jitter)
      U <- FZ(ZX)
    else {
      ## jitter observed values to avoid ties
      grain <- min(diff(sortunique(ZX)))/8
      jit <- runif(length(ZX), min=0, max=grain)
      sgn <- sample(c(-1L,1L), length(ZX), replace=TRUE)
      sgn[ZX==min(xxx)] <- 1L
      sgn[ZX==max(xxx)] <- -1L
      U <- FZ(ZX + sgn*jit)
    }

    ## Test uniformity
    result <- switch(test,
                     ks  = ks.test(U, "punif", ...),
                     cvm = cvm.test(U, "punif", ...),
                     ad = ad.test(U, "punif", ...))
    testname <- switch(test,
                       ks="Kolmogorov-Smirnov",
                       cvm="Cramer-Von Mises",
                       ad="Anderson-Darling")
    result$method <- paste("Spatial", testname, "test")
    result$data.name <-
      paste("predicted cdf of covariate", sQuote(paste(covname, collapse="")),
            "evaluated at data points of", sQuote(modelname))
    if(verbose)
      cat("done.\n")
    class(result) <- c("cdftest", class(result))
    attr(result, "prep") <- list(Zvalues = Zvalues, lambda = lambda,
                                 ZX = ZX, FZ = FZ, U = U, type = type)
    attr(result, "info") <- list(modelname = modelname, covname = covname)
    return(result)        
  }

  cdf.test.mppm

})

