##
## auc.R
##
##  Calculate ROC curve or area under it
##
## $Revision: 1.4 $ $Date: 2015/04/15 09:23:22 $

roc <- function(X, ...) { UseMethod("roc") }

roc.ppp <- function(X, covariate, ..., high=TRUE) {
  d <- spatialCDFframe(ppm(X), covariate, ...)
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
  d <- spatialCDFframe(ppm(Y), lambda, ...)
  U <- d$values$U
  ec <- ecdf(1-U) 
  p <- seq(0,1,length=1024)
  fobs <- ec(p)
  FZ <- d$values$FZ
  F1Z <- ewcdf(lambda[], lambda[]/sum(lambda))
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

auc <- function(X, ...) { UseMethod("auc") }

auc.ppp <- function(X, covariate, ..., high=TRUE) {
  d <- spatialCDFframe(ppm(X), covariate, ...)
  U <- d$values$U
  EU <- mean(U)
  result <- if(high) EU else (1 - EU) 
  return(result)
}

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

