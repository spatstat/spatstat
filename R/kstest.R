#
#  kstest.R
#
#  $Revision: 1.58 $  $Date: 2013/01/30 02:12:56 $
#
#

# --------- old -------------

ks.test.ppm <- function(...) {
  .Deprecated("kstest.ppm", package="spatstat")
  kstest.ppm(...)
}

# ---------------------------

kstest <- function(...) {
  UseMethod("kstest")
}

kstest.ppp <-
  function(X, covariate, ..., jitter=TRUE) {
    Xname <- short.deparse(substitute(X))
    covname <- singlestring(short.deparse(substitute(covariate)))
    if(is.character(covariate)) covname <- covariate
    if(!is.marked(X, dfok=TRUE)) {
      # unmarked
      model <- ppm(X)
      modelname <- "CSR"
    } else if(is.multitype(X)) {
      # multitype
      mf <- summary(X)$marks$frequency
      if(all(mf > 0)) {
        model <- ppm(X, ~marks)
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
    do.call("spatialCDFtest",
            resolve.defaults(list(model, covariate, test="ks"),
                             list(jitter=jitter),
                             list(...),
                             list(modelname=modelname,
                                  covname=covname, dataname=Xname)))
}

kstest.ppm <- 
  function(model, covariate, ..., jitter=TRUE) {
  modelname <- short.deparse(substitute(model))
  covname <- singlestring(short.deparse(substitute(covariate)))
  verifyclass(model, "ppm")
  if(is.character(covariate)) covname <- covariate
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call("spatialCDFtest",
          resolve.defaults(list(model, covariate, test="ks"),
                           list(jitter=jitter),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}

kstest.lpp <-
  function(X, covariate, ..., jitter=TRUE) {
    Xname <- short.deparse(substitute(X))
    covname <- singlestring(short.deparse(substitute(covariate)))
    if(is.character(covariate)) covname <- covariate
    if(!is.marked(X, dfok=TRUE)) {
      # unmarked
      model <- lppm(X)
      modelname <- "CSR"
    } else if(is.multitype(X)) {
      # multitype
      mf <- table(marks(X))
      if(all(mf > 0)) {
        model <- lppm(X, ~marks)
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
    do.call("spatialCDFtest",
            resolve.defaults(list(model, covariate, test="ks"),
                             list(jitter=jitter),
                             list(...),
                             list(modelname=modelname,
                                  covname=covname, dataname=Xname)))
}

kstest.lppm <- function(model, covariate, ..., jitter=TRUE) {
  modelname <- short.deparse(substitute(model))
  covname <- singlestring(short.deparse(substitute(covariate)))
  verifyclass(model, "lppm")
  if(is.character(covariate)) covname <- covariate
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call("spatialCDFtest",
          resolve.defaults(list(model, covariate, test="ks"),
                           list(jitter=jitter),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}


kstest.slrm <- function(model, covariate, ..., modelname=NULL, covname=NULL) {
  # get names
  if(is.null(modelname))
    modelname <- short.deparse(substitute(model))
  if(is.null(covname))
    covname <- short.deparse(substitute(covariate))
  dataname <- model$CallInfo$responsename
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
  FZX <- ewcdf(ZX)
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
  result <- ks.test(U, "punif", ...)

  # modify the 'htest' entries
  result$method <- paste("Spatial Kolmogorov-Smirnov test of",
                         "inhomogeneous Poisson process",
                         "in two dimensions")
  result$data.name <-
    paste("covariate", sQuote(paste(covname, collapse="")),
          "evaluated at points of", sQuote(dataname), "\n\t",
          "and transformed to uniform distribution under model",
          sQuote(modelname))
  # additional class 'kstest'
  class(result) <- c("kstest", class(result))
  attr(result, "prep") <-
    list(Zvalues=covvalu, ZX=ZX, FZ=FZ, FZX=ecdf(ZX), U=U)
  attr(result, "info") <- list(modelname=modelname, covname=covname,
                               dataname=dataname, csr=FALSE)
  return(result)        
}

#.............  helper functions ........................#

spatialCDFtest <- function(model, covariate, test, ...,
                           dimyx=NULL, eps=NULL,
                           jitter=TRUE, 
                           modelname=NULL, covname=NULL, dataname=NULL) {
  if(!is.poisson(model))
    stop("Only implemented for Poisson point process models")
  # conduct test based on comparison of CDF's of covariate values
  test <- pickoption("test", test, c(ks="ks"))
  # compute the essential data
  fra <- spatialCDFframe(model, covariate,
                         dimyx=dimyx, eps=eps,
                         jitter=jitter, modelname=modelname,
                         covname=covname, dataname=dataname)
  values <- fra$values
  info   <- fra$info
  # Test uniformity of transformed values
  U <- values$U
  switch(test,
         ks={ result <- ks.test(U, "punif", ...) },
         stop("Unrecognised test option"))

  # modify the 'htest' entries
  csr <- info$csr
  result$method <- paste("Spatial Kolmogorov-Smirnov test of",
                         if(csr) "CSR" else "inhomogeneous Poisson process",
                         "in", info$spacename)
  result$data.name <-
    paste("covariate", sQuote(singlestring(info$covname)),
          "evaluated at points of", sQuote(info$dataname), "\n\t",
          "and transformed to uniform distribution under",
          if(csr) info$modelname else sQuote(info$modelname))
  
  # additional class 'kstest'
  class(result) <- c("kstest", class(result))
  attr(result, "frame") <- fra
  return(result)        
}

spatialCDFframe <- function(model, covariate, ...) {
  # evaluate CDF of covariate values at data points and at pixels
  stuff <- evalCovar(model, covariate, ...)
  # extract 
  values <- stuff$values
  info   <- stuff$info
  Zvalues <- values$Zvalues
  lambda  <- values$lambda
  weights <- values$weights
  ZX      <- values$ZX
  # compute empirical cdf of Z values at points of X
  FZX <- ecdf(ZX)
  # form weighted cdf of Z values in window
  wts <- lambda * weights
  FZ <- ewcdf(Zvalues, wts/sum(wts))
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

  # pack up
  stuff$values$FZ  <- FZ
  stuff$values$FZX <- FZX
  stuff$values$U   <- U
  class(stuff) <- "spatialCDFframe"
  return(stuff)
}

plot.kstest <- function(x, ..., style=c("cdf", "PP", "QQ"),
                        lwd=par("lwd"), col=par("col"), lty=par("lty"),
                        lwd0=lwd, col0=col, lty0=lty) {
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
           do.call("plot.default",
                   resolve.defaults(
                                    list(x=qZ, y=pZ, type="l"),
                                    list(...),
                                    list(lwd=lwd0, col=col0, lty=lty0),
                                    list(xlab=info$covname, ylab="probability",
                                         main=main)))
           plot(FZX, add=TRUE, do.points=FALSE, lwd=lwd, col=col, lty=lty)
         },
         PP={
           # plot FZX o (FZ)^{-1}
           pX <- get("y", environment(FZX))
           qX <- get("x", environment(FZX))
           p0 <- FZ(qX)
           do.call("plot.default",
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
           do.call("plot.default",
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
