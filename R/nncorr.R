#
# nncorr.R
#
# $Revision: 1.12 $  $Date: 2019/01/22 03:08:57 $
#

nnmean <- function(X, k=1, na.action="warn") {
  stopifnot(is.ppp(X))
  if(!is.marked(X, na.action=na.action))
    stop("X must be a marked point pattern", call.=FALSE)
  if(k %% 1 != 0 || length(k) != 1 || k <= 0)
    stop("k should be a single integer greater than 0", call.=FALSE)
  m <- numeric.columns(marks(X), logical=TRUE, others="na")
  ## default result
  nana <- rep(NA_real_, ncol(m))
  ans <- rbind(unnormalised=nana,
               normalised=nana)
  ## 
  if(all(is.na(m))) {
    warning("non-numeric marks; results are NA", call.=FALSE)
  } else if(k >= npoints(X)) {
    warning(paste("Not enough points to compute k-th nearest neighbours",
                  paste0(paren(paste0("n = ", npoints(X), ", k = ", k)), ";"),
                  "results are NA"),
            call.=FALSE)
  } else {
    nnid <- nnwhich(X, k=k)
    ok <- (nndist(X, k=k) <= bdist.points(X))
    if(!any(ok, na.rm=TRUE)) {
      warning("insufficient data remaining after border correction; results are NA")
    } else {
      numer <- sapply(as.data.frame(m[nnid[ok], ]), mean, na.rm=TRUE)
      denom <- sapply(as.data.frame(m),             mean, na.rm=TRUE)
      ans <- rbind(unnormalised=numer,
                   normalised  =numer/denom)
    }
  }
  if(ncol(ans) == 1) ans <- ans[,1,drop=TRUE]
  return(ans)
}

nnvario <- local({

  nnvario <- function(X, k=1, na.action="warn") {
    stopifnot(is.ppp(X))
    if(!is.marked(X, na.action=na.action))
      stop("X must be a marked point pattern", call.=FALSE)
    m <- numeric.columns(marks(X), logical=TRUE, others="na")
    if(all(is.na(m))) warning("non-numeric marks; results are NA", call.=FALSE)
    ans <- nncorr(X %mark% m, sqdif, k=k, denominator=diag(var(m)),
                  na.action="ignore")
    return(ans)
  }
  sqdif <- function(m1,m2) { ((m1-m2)^2)/2 }

  nnvario
})


nncorr <- function(X, f = function(m1,m2) { m1 * m2},
                   k=1,
                   ...,
                   use = "all.obs",
                   method = c("pearson", "kendall", "spearman"),
                   denominator=NULL,
                   na.action="warn") {
  stopifnot(is.ppp(X))
  if(!is.marked(X, na.action=na.action))
    stop("X must be a marked point pattern", call.=FALSE)
  if(k %% 1 != 0 || length(k) != 1 || k <= 0)
    stop("k should be a single integer greater than 0", call.=FALSE)
  if(k >= npoints(X))
    stop("Not enough points to compute k-th nearest neighbours")
  
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
    if(is.function(f)) f <- rep.int(list(f), nv)
    if(length(denominator) <= 1) denominator <- rep.int(list(denominator), nv)
    #
    result <- matrix(NA, nrow=3, ncol=nv)
    outnames <- c("unnormalised", "normalised", "correlation")
    dimnames(result) <- list(outnames, colnames(m))
    for(j in 1:nv) {
      mj <- m[,j, drop=FALSE]
      denj <- denominator[[j]]
      nncj <- nncorr(X %mark% mj, f=f[[j]], k=k, use=use, method=method,
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
  nn <- nnwhich(X, k=k)
  ok <- (nndist(X, k=k) <= bdist.points(X))
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
  
