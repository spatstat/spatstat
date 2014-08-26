#
# cdf.test.mppm.R
#
# $Revision: 1.12 $  $Date: 2014/06/10 02:47:08 $
#
cdf.test.mppm <- function(model, covariate,
                          test=c("ks", "cvm", "ad"), ..., verbose=TRUE,
                         interpolate=FALSE, fast=TRUE, jitter=TRUE) {
  modelname <- short.deparse(substitute(model))
  covname <- short.deparse(substitute(covariate))
  test <- match.arg(test)
  stopifnot(is.mppm(model))
  if(!is.poisson.mppm(model))
    stop("Only implemented for Poisson models")
  pixelvalues <- function(z) { as.vector(as.matrix(z)) }
  # extract things from model
  data  <- model$data
  npat  <- model$npat
  Y     <- data.mppm(model)
  if(fast) {
    # extract original quadrature schemes and convert to point patterns
    QQ  <- quad.mppm(model)
    PP  <- lapply(QQ, union.quad)
    Zweights <- lapply(QQ, w.quad)
  } else
    Zweights <- list()
  # `evaluate' covariate
  if(verbose)
    cat("Extracting covariate...")
  if(is.character(covariate)) {
    # extract covariate with this name from data used to fit model
    if(!(covariate %in% names(data)))
      stop(paste("Model does not contain a covariate called",
                 dQuote(covariate)))
    covname <- covariate
    covariate <- data[, covname, drop=TRUE]
  } else if(inherits(covariate, "listof")) {
    if(length(covariate) != npat)
      stop(paste("Length of list of covariate values does not match",
                 "number of point patterns in data of original model"))
  } else if(is.hyperframe(covariate)) {
    # extract first column
    covariate <- covariate[,1, drop=TRUE]
    if(length(covariate) != npat)
      stop(paste("Number of rows of covariate hyperframe does not match",
                 "number of point patterns in data of original model"))
  } else if(is.function(covariate) || is.im(covariate)) {
    # replicate to make a list
    covariate <- rep(list(covariate), npat)
    class(covariate) <- c("listof", class(covariate))
  } else     
  stop(paste("Format of argument", sQuote("covariates"), "not understood"))
  if(verbose)
    cat("done.\nComputing statistics for each pattern...")

  # compile information for test from each row
  Zvalues <- ZX <- Win <- list()
  for(i in 1:npat) {
    if(verbose) progressreport(i, npat)
    XI <- Y[[i]]
    if(fast)
      PI <- PP[[i]]
    else
      WI <- XI$window
    covariateI <- covariate[[i]]
    if(is.im(covariateI)) {
      type <- "im"
      # evaluate at data points
      ZXI <-
        if(interpolate) interp.im(covariateI, XI$x, XI$y)
        else covariateI[XI]
      if(fast) {
        # covariate values for quadrature points
        ZI <- covariateI[PI]
      } else {
        # covariate image inside window
        ZI <- covariateI[WI, drop=FALSE]
        # corresponding mask
        WI <- as.owin(ZI)
        # pixel areas 
        Zweights[[i]] <- rep(WI$xstep * WI$ystep, prod(WI$dim))
      }
    } else if(is.function(covariateI)) {
      type <- "function"
      # evaluate exactly at data points
      ZXI <- covariateI(XI$x, XI$y)
      if(fast) {
        # covariate values for quadrature points
        ZI <- covariateI(PI$x, PI$y)
      } else {
        # window
        WI <- as.mask(WI)
        # covariate image inside window
        ZI <- as.im(covariateI, W=WI)
        # pixel areas 
        Zweights[[i]] <- rep(WI$xstep * WI$ystep, prod(WI$dim))
      }
    } else
    stop("covariate should be an image or a function(x,y)")
    ZX[[i]] <- ZXI
    if(fast)
      Zvalues[[i]] <- ZI      
    else {
      Win[[i]] <- WI
      # values of covariate in window
      Zvalues[[i]] <- pixelvalues(ZI)
    }
  }

  if(verbose)
    cat("done.\nComputing predicted intensity...")

  # compute predicted intensities
  trend <-
    if(fast)
      fitted(model, type="trend")
    else
      predict(model, type="trend", locations=Win, verbose=verbose)$trend
  
  if(verbose)
    cat("done.\nExtracting...")
  # extract relevant values
  lambda <- if(fast) trend else lapply(trend, pixelvalues)
  if(verbose)
    cat("done.\nPerforming test...")
  
  # flatten to vectors
  lambda <- unlist(lambda)
  Zweights <- unlist(Zweights)
  Zvalues <- unlist(Zvalues)
  ZX      <- unlist(ZX)
  if(length(lambda) != length(Zvalues))
    stop("Internal error: mismatch between predicted values and Z values")
  if(length(Zvalues) != length(Zweights))
    stop("Internal error: mismatch between Z values and Z weights")
  lambda <- lambda * Zweights
  
  # form weighted cdf of Z values in window
  FZ <- ewcdf(Zvalues, lambda/sum(lambda))
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
  # evaluate at data points
  if(!jitter)
    U <- FZ(ZX)
  else {
    # jitter observed values to avoid ties
    grain <- min(diff(sort(unique(ZX))))/8
    jit <- runif(length(ZX), min=0, max=grain)
    sgn <- sample(c(-1,1), length(ZX), replace=TRUE)
    sgn[ZX==min(xxx)] <- 1
    sgn[ZX==max(xxx)] <- -1
    U <- FZ(ZX + sgn*jit)
  }

  # Test uniformity
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
