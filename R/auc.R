##
## auc.R
##
##  Calculate ROC curve or area under it
##
## $Revision: 1.6 $ $Date: 2016/11/10 01:08:04 $

roc <- function(X, ...) { UseMethod("roc") }

roc.ppp <- function(X, covariate, ..., high=TRUE) {
  nullmodel <- ppm(X)
  result <- rocData(covariate, nullmodel, ..., high=high)
  return(result)
}

roc.lpp <- function(X, covariate, ..., high=TRUE) {
  nullmodel <- lppm(X)
  result <- rocData(covariate, nullmodel, ..., high=high)
  return(result)
}

rocData <- function(covariate, nullmodel, ..., high=TRUE) {
  d <- spatialCDFframe(nullmodel, covariate, ...)
  U <- d$values$U
  ec <- if(high) ecdf(1-U) else ecdf(U)
  p <- seq(0,1,length=1024)
  df <- data.frame(p=p, fobs=ec(p), fnull=p)
  result <- fv(df,
               argu="p",
               ylab=quote(roc(p)),
               valu="fobs",
               desc=c("fraction of area",
                      "observed fraction of points",
                      "expected fraction if no effect"),
               fname="roc")
  fvnames(result, ".") <- c("fobs", "fnull")
  return(result)
}

roc.ppm <- function(X, ...) {
  stopifnot(is.ppm(X))
  model <- X
  lambda <- predict(model, ...)
  Y <- data.ppm(model)
  nullmodel <- ppm(Y)
  result <- rocModel(lambda, nullmodel, ...)
  return(result)
}

roc.kppm <- function(X, ...) {
  stopifnot(is.kppm(X))
  model <- as.ppm(X)
  lambda <- predict(model, ...)
  Y <- data.ppm(model)
  nullmodel <- ppm(Y)
  result <- rocModel(lambda, nullmodel, ...)
  return(result)
}

roc.lppm <- function(X, ...) {
  stopifnot(is.lppm(X))
  model <- X
  lambda <- predict(model, ...)
  Y <- X$X
  nullmodel <- lppm(Y)
  result <- rocModel(lambda, nullmodel, ...)
  return(result)
}

rocModel <- function(lambda, nullmodel, ..., high) {
  if(!missing(high))
    warning("Argument 'high' is ignored when computing ROC for a fitted model")
  d<- spatialCDFframe(nullmodel, lambda, ...) 
  U <- d$values$U
  ec <- ecdf(1-U) 
  p <- seq(0,1,length=1024)
  fobs <- ec(p)
  FZ <- d$values$FZ
  lambdavalues <- if(is.im(lambda)) lambda[] else unlist(lapply(lambda, "["))
  F1Z <- ewcdf(lambdavalues, lambdavalues/sum(lambdavalues))    
  pZ <- get("y", environment(FZ))
  qZ <- get("x", environment(FZ))
  FZinverse <- approxfun(pZ, qZ, rule=2)
  ftheo <- 1 - F1Z(FZinverse(1-p))
  df <- data.frame(p=p, fobs=fobs, ftheo=ftheo, fnull=p)
  result <- fv(df,
               argu="p",
               ylab=quote(roc(p)),
               valu="fobs",
               fmla = . ~ p,
               desc=c("fraction of area",
                 "observed fraction of points",
                 "expected fraction of points",
                 "expected fraction if no effect"),
               fname="roc")
  fvnames(result, ".") <- c("fobs", "ftheo", "fnull")
  return(result)
}

#    ......................................................

auc <- function(X, ...) { UseMethod("auc") }

auc.ppp <- function(X, covariate, ..., high=TRUE) {
  d <- spatialCDFframe(ppm(X), covariate, ...)
  U <- d$values$U
  EU <- mean(U)
  result <- if(high) EU else (1 - EU) 
  return(result)
}

auc.lpp <- function(X, covariate, ..., high=TRUE) {
  d <- spatialCDFframe(lppm(X), covariate, ...)
  U <- d$values$U
  EU <- mean(U)
  result <- if(high) EU else (1 - EU) 
  return(result)
}

auc.kppm <- function(X, ...) { auc(as.ppm(X), ...) }

auc.ppm <- function(X, ...) {
  model <- X
  if(is.multitype(model)) {
    # cheat
    ro <- roc(model, ...)
    aobs <- with(ro, mean(fobs))
    atheo <- with(ro, mean(ftheo))
  } else if(is.stationary(model)) {
    aobs <- atheo <- 1/2
  } else {
    lambda <- intensity(model)
    Fl <- ecdf(lambda[])
    lambda <- as.im(lambda, Window(model))
    X <- data.ppm(model)
    lamX <- lambda[X]
    aobs <- mean(Fl(lamX))
    atheo <- mean(lambda[] * Fl(lambda[]))/mean(lambda)
  }
  result <- c(aobs, atheo)
  names(result) <- c("obs", "theo")
  return(result)
}

auc.lppm <- function(X, ...) {
  stopifnot(inherits(X, "lppm"))
  model <- X
  if(is.multitype(model)) {
    # cheat
    ro <- roc(model, ...)
    aobs <- with(ro, mean(fobs))
    atheo <- with(ro, mean(ftheo))
  } else {
    lambda <- predict(model, ...)
    Fl <- ecdf(lambda[])
    lamX <- lambda[model$X]
    aobs <- mean(Fl(lamX))
    atheo <- mean(lambda[] * Fl(lambda[]))/mean(lambda)
  }
  result <- c(aobs, atheo)
  names(result) <- c("obs", "theo")
  return(result)
}

