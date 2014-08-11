#
# nncorr.R
#
# $Revision: 1.7 $  $Date: 2011/10/11 10:45:24 $
#

nnmean <- function(X) {
  stopifnot(is.ppp(X) && is.marked(X))
  m <- numeric.columns(marks(X), logical=TRUE, others="na")
  nv <- ncol(m)
  nnid <- nnwhich(X)
  ok <- (nndist(X) <= bdist.points(X))
  if(!any(ok))
    stop("Insufficient data")
  numer <- unlist(lapply(as.data.frame(m[nnid[ok], ]), mean, na.rm=TRUE))
  denom <- unlist(lapply(as.data.frame(m),             mean, na.rm=TRUE))
  ans <- rbind(unnormalised=numer,
               normalised  =numer/denom)
  return(ans)
}

nnvario <- function(X) {
  stopifnot(is.ppp(X) && is.marked(X))
  m <- numeric.columns(marks(X), logical=TRUE, others="na")
  f <- function(m1,m2) { ((m1-m2)^2)/2 }
  ans <- nncorr(X %mark% m, f, denominator=diag(var(m)))
  return(ans)
}

nncorr <- function(X, f = function(m1,m2) { m1 * m2}, ...,
                   use = "all.obs",
                   method = c("pearson", "kendall", "spearman"),
                   denominator=NULL) {
  stopifnot(is.ppp(X) && is.marked(X))
  m <- as.data.frame(marks(X))
  nv <- ncol(m)
  if(nv == 1) colnames(m) <- ""
  #
  if(missing(method) || is.null(method))
    method <- "pearson"
  # 
  if(missing(f)) f <- NULL
  if(!is.null(f) && !is.function(f)) {
    if(nv == 1) stop("f should be a function")
    # could be a list of functions
    if(!(is.list(f) && all(unlist(lapply(f, is.function)))))
      stop("f should be a function or a list of functions")
    if(length(f) != nv)
      stop("Length of list f does not match number of mark variables")
  }
  # optional denominator(s)
  if(!is.null(denominator) && !(length(denominator) %in% c(1, nv)))
    stop("Denominator has incorrect length")
  # multi-dimensional case
  if(nv > 1) {
    # replicate things
    if(is.function(f)) f <- rep(list(f), nv)
    if(length(denominator) <= 1) denominator <- rep(list(denominator), nv)
    #
    result <- matrix(NA, nrow=3, ncol=nv)
    outnames <- c("unnormalised", "normalised", "correlation")
    dimnames(result) <- list(outnames, colnames(m))
    for(j in 1:nv) {
      mj <- m[,j, drop=FALSE]
      denj <- denominator[[j]]
      nncj <- nncorr(X %mark% mj, f=f[[j]], use=use, method=method,
                     denominator=denj)
      kj <- length(nncj)
      result[1:kj,j] <- nncj
    }
    if(all(is.na(result[3, ]))) result <- result[1:2, ]
    return(result)
  }
  # one-dimensional
  m <- m[,1,drop=TRUE]
  # select 'f' appropriately for X
  chk <- check.testfun(f, X=X)
  f     <- chk$f
  ftype <- chk$ftype
  # denominator
  Efmm <-
    if(!is.null(denominator)) denominator else 
    switch(ftype,
           mul={ 
             mean(m)^2
           },
           equ={
             sum(table(m)^2)/length(m)^2
           },
           general={
             mean(outer(m, m, f, ...))
           })
  # border method
  nn <- nnwhich(X)
  ok <- (nndist(X) <= bdist.points(X))
  if(!any(ok))
    stop("Insufficient data")
  mY <- m[nn[ok]]
  mX <- m[ok]
  Efmk <- switch(ftype,
                 mul = {
                   mean(mX * mY, ...)
                 },
                 equ = {
                   mean(mX == mY, ...)
                 }, 
                 general = {
                   mean(f(mX, mY, ...))
                 })
  #
  answer <- c(unnormalised=Efmk,
              normalised=Efmk/Efmm)
  if(ftype == "mul") {
    classic <- cor(mX, mY, use=use, method=method)
    answer <- c(answer, correlation=classic)
  }
  return(answer)
}
  
