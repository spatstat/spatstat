#
#  cdftest.R
#
#  $Revision: 2.19 $  $Date: 2018/04/19 02:00:41 $
#
#

cdf.test <- function(...) {
  UseMethod("cdf.test")
}

cdf.test.ppp <-
  function(X, covariate, test=c("ks", "cvm", "ad"), ...,
           interpolate=TRUE, jitter=TRUE) {
    Xname <- short.deparse(substitute(X))
    covname <- singlestring(short.deparse(substitute(covariate)))
    test <- match.arg(test)
    if(is.character(covariate)) covname <- covariate
    if(!is.marked(X, dfok=TRUE)) {
      # unmarked
      model <- ppm(X)
      modelname <- "CSR"
    } else if(is.multitype(X)) {
      # multitype
      mf <- summary(X)$marks$frequency
      if(all(mf > 0)) {
        model <- ppm(X ~marks)
        modelname <- "CSRI"
      } else {
        warning("Ignoring marks, because some mark values have zero frequency")
        X <- unmark(X)
        model <- ppm(X)
        modelname <- "CSR"
      } 
    } else {
      # marked - general case
      X <- unmark(X)
      warning("marks ignored")
      model <- ppm(X)
      modelname <- "CSR"
    }
    do.call(spatialCDFtest,
            resolve.defaults(list(model, covariate, test=test),
                             list(interpolate=interpolate, jitter=jitter),
                             list(...),
                             list(modelname=modelname,
                                  covname=covname, dataname=Xname)))
}

cdf.test.ppm <- 
  function(model, covariate, test=c("ks", "cvm", "ad"), ...,
           interpolate=TRUE, jitter=TRUE, nsim=99, verbose=TRUE) {
  modelname <- short.deparse(substitute(model))
  covname <- singlestring(short.deparse(substitute(covariate)))
  test <- match.arg(test)
  verifyclass(model, "ppm")
  if(is.character(covariate)) covname <- covariate
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call(spatialCDFtest,
          resolve.defaults(list(model, covariate, test=test),
                           list(interpolate=interpolate, jitter=jitter,
                                nsim=nsim, verbose=verbose),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}

cdf.test.lpp <-
  function(X, covariate, test=c("ks", "cvm", "ad"), ...,
           interpolate=TRUE, jitter=TRUE) {
    Xname <- short.deparse(substitute(X))
    covname <- singlestring(short.deparse(substitute(covariate)))
    test <- match.arg(test)
    if(is.character(covariate)) covname <- covariate
    if(!is.marked(X, dfok=TRUE)) {
      # unmarked
      model <- lppm(X)
      modelname <- "CSR"
    } else if(is.multitype(X)) {
      # multitype
      mf <- table(marks(X))
      if(all(mf > 0)) {
        model <- lppm(X ~ marks)
        modelname <- "CSRI"
      } else {
        warning("Ignoring marks, because some mark values have zero frequency")
        X <- unmark(X)
        model <- ppm(X)
        modelname <- "CSR"
      } 
    } else {
      # marked - general case
      X <- unmark(X)
      warning("marks ignored")
      model <- ppm(X)
      modelname <- "CSR"
    }
    do.call(spatialCDFtest,
            resolve.defaults(list(model, covariate, test=test),
                             list(interpolate=interpolate, jitter=jitter),
                             list(...),
                             list(modelname=modelname,
                                  covname=covname, dataname=Xname)))
}

cdf.test.lppm <- function(model, covariate,
                          test=c("ks", "cvm", "ad"),
                          ..., interpolate=TRUE, jitter=TRUE,
                          nsim=99, verbose=TRUE) {
  modelname <- short.deparse(substitute(model))
  covname <- singlestring(short.deparse(substitute(covariate)))
  test <- match.arg(test)
  verifyclass(model, "lppm")
  if(is.character(covariate)) covname <- covariate
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call(spatialCDFtest,
          resolve.defaults(list(model, covariate, test=test),
                           list(interpolate=interpolate, jitter=jitter,
                                nsim=nsim, verbose=verbose),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}


cdf.test.slrm <- function(model, covariate,
                          test=c("ks", "cvm", "ad"), ...,
                          modelname=NULL, covname=NULL) {
  # get names
  if(is.null(modelname))
    modelname <- short.deparse(substitute(model))
  if(is.null(covname))
    covname <- short.deparse(substitute(covariate))
  dataname <- model$CallInfo$responsename
  test <- match.arg(test)
  #
  stopifnot(is.slrm(model))
  stopifnot(is.im(covariate))
  # extract data
  prob <- fitted(model)
  covim <- as.im(covariate, W=as.owin(prob))
  probvalu <- as.matrix(prob)
  covvalu  <- as.matrix(covim)
  ok <- !is.na(probvalu) & !is.na(covvalu)
  probvalu <- as.vector(probvalu[ok])
  covvalu <- as.vector(covvalu[ok])
  # compile weighted cdf's
  FZ <- ewcdf(covvalu, probvalu/sum(probvalu))
  X <- model$Data$response
  ZX <- safelookup(covim, X)
  # Ensure support of cdf includes the range of the data
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
  # make piecewise linear approximation of cdf
  FZ <- approxfun(xxx, yyy, rule=2)
  # now apply cdf
  U <- FZ(ZX)
  # Test uniformity of transformed values
  result <- switch(test,
                   ks  = ks.test(U, "punif", ...),
                   cvm = cvm.test(U, "punif", ...),
                   ad = ad.test(U, "punif", ...))
  testname <- switch(test,
                     ks="Kolmogorov-Smirnov",
                     cvm="Cramer-Von Mises",
                     ad="Anderson-Darling")

  # modify the 'htest' entries
  result$method <- paste("Spatial", testname, "test of",
                         "inhomogeneous Poisson process",
                         "in two dimensions")
  result$data.name <-
    paste("covariate", sQuote(paste(covname, collapse="")),
          "evaluated at points of", sQuote(dataname),
          "\n     and transformed to uniform distribution under",
          sQuote(modelname))
  # additional class 'cdftest'
  class(result) <- c("cdftest", class(result))
  attr(result, "prep") <-
    list(Zvalues=covvalu, ZX=ZX, FZ=FZ, FZX=ecdf(ZX), U=U)
  attr(result, "info") <- list(modelname=modelname, covname=covname,
                               dataname=dataname, csr=FALSE)
  return(result)        
}

#.............  helper functions ........................#

spatialCDFtest <- function(model, covariate, test=c("ks", "cvm", "ad"),
                           ...,
                           dimyx=NULL, eps=NULL,
                           interpolate=TRUE, jitter=TRUE, nsim=99, verbose=TRUE,
                           modelname=NULL, covname=NULL, dataname=NULL) {
  ## conduct test based on comparison of CDF's of covariate values
  test <- match.arg(test)
  ## compute the essential data
  fra <- spatialCDFframe(model, covariate,
                         dimyx=dimyx, eps=eps,
                         interpolate=interpolate, jitter=jitter,
                         modelname=modelname,
                         covname=covname, dataname=dataname)
  ## calculate the test statistic
  result <- spatialCDFtestCalc(fra, test=test, ...)

  if(is.poisson(model))
    return(result)

  ## Gibbs model: perform Monte Carlo test
  result$poisson.p.value <- pobs <- result$p.value
  result$poisson.statistic <- tobs <- result$statistic
  Xsim <- simulate(model, nsim=nsim, progress=verbose)
  sim.pvals <- sim.stats <- numeric(nsim)
  if(verbose) {
    cat("Processing.. ")
    state <- list()
  }
  for(i in seq_len(nsim)) {
    model.i <- update(model, Xsim[[i]])
    fra.i <- spatialCDFframe(model.i, covariate,
                             dimyx=dimyx, eps=eps,
                             interpolate=interpolate, jitter=jitter,
                             modelname=modelname,
                             covname=covname, dataname=dataname)
    res.i <- spatialCDFtestCalc(fra.i, test=test, ..., details=FALSE)
    sim.pvals[i] <- res.i$p.value
    sim.stats[i] <- res.i$statistic
    if(verbose) state <- progressreport(i, nsim, state=state)
  }
  if(verbose) cat("Done.\n")
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
  ## modify the 'htest' entries
  testname <- switch(test,
                     ks="Kolmogorov-Smirnov",
                     cvm="Cramer-Von Mises",
                     ad="Anderson-Darling")
  result$method <-
    paste("Monte Carlo spatial", testname, "test",
          "of Gibbs process in", fra$info$spacename)
  return(result)        
}

spatialCDFtestCalc <- function(fra, test=c("ks", "cvm", "ad"), ...,
                               details=TRUE) {
  test <- match.arg(test)
  values <- fra$values
  info   <- fra$info
  ## Test uniformity of transformed values
  U <- values$U
  result <- switch(test,
                   ks  = ks.test(U, "punif", ...),
                   cvm = cvm.test(U, "punif", ...),
                   ad = ad.test(U, "punif", ...))

  # shortcut for internal use only
  if(!details) 
    return(result)
  
  ## add a full explanation, internal data, etc.
  
  ## modify the 'htest' entries
  csr    <- info$csr
  ispois <- info$ispois
  modelname <-
    if(csr) "CSR" else
    if(ispois) "inhomogeneous Poisson process" else "Gibbs process"
  testname <- switch(test,
                     ks="Kolmogorov-Smirnov",
                     cvm="Cramer-Von Mises",
                     ad="Anderson-Darling")
  result$method <-
    paste("Spatial", testname, "test of", modelname, "in", info$spacename)
  result$data.name <-
    paste("covariate", sQuote(singlestring(info$covname)),
          "evaluated at points of", sQuote(info$dataname), 
          "\n     and transformed to uniform distribution under",
          if(csr) info$modelname else sQuote(info$modelname))

  ## include internal data
  attr(result, "frame") <- fra

  ## additional class 'cdftest'
  class(result) <- c("cdftest", class(result))
  return(result)        
}

spatialCDFframe <- function(model, covariate, ..., jitter=TRUE) {
  # evaluate CDF of covariate values at data points and at pixels
  stuff <- evalCovar(model, covariate, ..., jitter=jitter)
  # extract 
  values <- stuff$values
#  info   <- stuff$info
  Zvalues <- values$Zvalues
  lambda  <- values$lambda
  weights <- values$weights
  ZX      <- values$ZX
  # compute empirical cdf of Z values at points of X
  FZX <- ecdf(ZX)
  # form weighted cdf of Z values in window
  wts <- lambda * weights
  sumwts <- sum(wts)
  FZ <- ewcdf(Zvalues, wts/sumwts)
  # Ensure support of cdf includes the range of the data
  xxx <- knots(FZ)
  yyy <- FZ(xxx)
  minZX <- min(ZX, na.rm=TRUE)
  minxxx <- min(xxx, na.rm=TRUE)
  if(minxxx > minZX) {
    xxx <- c(minZX, xxx)
    yyy <- c(0, yyy)
  }
  maxZX <- max(ZX, na.rm=TRUE)
  maxxxx <- max(xxx, na.rm=TRUE)
  if(maxxxx < maxZX) {
    xxx <- c(xxx, maxZX)
    yyy <- c(yyy, 1)
  }
  # make piecewise linear approximation of cdf
  FZ <- approxfun(xxx, yyy, rule=2)
  # now apply cdf
  U <- FZ(ZX)

  if(jitter) {
    ## Z values have already been jittered, but this does not guarantee
    ## that U values are distinct
    nU <- length(U)
    U <- U + runif(nU, -1, 1)/max(100, 2*nU)
    U <- pmax(0, pmin(1, U))
  }
  
  # pack up
  stuff$values$FZ  <- FZ
  stuff$values$FZX <- FZX
  stuff$values$U   <- U
  stuff$values$EN <- sumwts  ## integral of intensity = expected number of pts
  class(stuff) <- "spatialCDFframe"
  return(stuff)
}

plot.cdftest <- function(x, ..., style=c("cdf", "PP", "QQ"),
                        lwd=par("lwd"), col=par("col"), lty=par("lty"),
                        lwd0=lwd, col0=2, lty0=2,
                        do.legend=TRUE) {
  style <- match.arg(style)
  fram <- attr(x, "frame")
  if(!is.null(fram)) {
    values <- fram$values
    info <- fram$info
  } else {
    # old style
    values <- attr(x, "prep")
    info <- attr(x, "info")
  }
  # cdf of covariate Z over window 
  FZ <- values$FZ
  # cdf of covariate values at data points
  FZX <- values$FZX
  # blurb
  covname <- info$covname
  covdescrip <- switch(covname,
                       x="x coordinate",
                       y="y coordinate",
                       paste("covariate", dQuote(covname)))
  # plot it
  switch(style,
         cdf={
           # plot both cdf's superimposed
           qZ <- get("x", environment(FZ))
           pZ <- get("y", environment(FZ))
           main <- c(x$method,
                     paste("based on distribution of", covdescrip),
                     paste("p-value=", signif(x$p.value, 4)))
           do.call(plot.default,
                   resolve.defaults(
                                    list(x=qZ, y=pZ, type="l"),
                                    list(...),
                                    list(lwd=lwd0, col=col0, lty=lty0),
                                    list(xlab=info$covname, ylab="probability",
                                         main=main)))
           plot(FZX, add=TRUE, do.points=FALSE, lwd=lwd, col=col, lty=lty)
           if(do.legend) 
             legend("topleft", c("observed", "expected"),
                    lwd=c(lwd,lwd0),
                    col=c(col2hex(col), col2hex(col0)),
                    lty=c(lty2char(lty),lty2char(lty0)))
         },
         PP={
           # plot FZX o (FZ)^{-1}
           pX <- get("y", environment(FZX))
           qX <- get("x", environment(FZX))
           p0 <- FZ(qX)
           do.call(plot.default,
                   resolve.defaults(
                                    list(x=p0, y=pX),
                                    list(...),
                                    list(col=col),
                                    list(xlim=c(0,1),
                                         ylim=c(0,1),
                                         xlab="Theoretical probability",
                                         ylab="Observed probability",
                                         main="")))
           abline(0,1, lwd=lwd0, col=col0, lty=lty0)           
         },
         QQ={
           # plot (FZX)^{-1} o FZ
           pZ <- get("y", environment(FZ))
           qZ <- get("x", environment(FZ))
           FZinverse <- approxfun(pZ, qZ, rule=2)
           pX <- get("y", environment(FZX))
           qX <- get("x", environment(FZX))
           qZX <- FZinverse(pX)
           Zrange <- range(qZ, qX, qZX)
           xlab <- paste("Theoretical quantile of", covname)
           ylab <- paste("Observed quantile of", covname)
           do.call(plot.default,
                   resolve.defaults(
                                    list(x=qZX, y=qX),
                                    list(...),
                                    list(col=col),
                                    list(xlim=Zrange, ylim=Zrange,
                                         xlab=xlab, ylab=ylab,
                                         main="")))
           abline(0,1, lwd=lwd0, col=col0, lty=lty0)           
         })
  return(invisible(NULL))
}
