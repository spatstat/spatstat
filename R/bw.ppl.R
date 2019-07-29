#'
#'   bw.ppl.R
#'
#'   Likelihood cross-validation for kernel smoother of point pattern
#'
#'   bw.ppl    class ppp
#'   bw.lppl   class lpp
#' 
#'   $Revision: 1.10 $ $Date: 2019/07/29 09:25:20 $
#'

bw.ppl <- function(X, ..., srange=NULL, ns=16, sigma=NULL,
                   weights=NULL, shortcut=FALSE) {
  stopifnot(is.ppp(X))
  if(!is.null(sigma)) {
    stopifnot(is.numeric(sigma) && is.vector(sigma))
    ns <- length(sigma)
  } else {
    if(!is.null(srange)) check.range(srange) else {
      nnd <- nndist(X)
      srange <- c(min(nnd[nnd > 0]), diameter(as.owin(X))/2)
    }
    sigma <- geomseq(from=srange[1L], to=srange[2L], length.out=ns)
  }
  cv <- numeric(ns)
  if(shortcut) {
    for(i in 1:ns) {
      si <- sigma[i]
      lamx <- density(X, sigma=si, at="points", leaveoneout=TRUE,
                    weights=weights, ...)
      cv[i] <- sum(log(lamx))
    }
  } else {
    for(i in 1:ns) {
      si <- sigma[i]
      lamx <- density(X, sigma=si, at="points", leaveoneout=TRUE,
                      weights=weights, ...)
      lam <- density(X, sigma=si,
                     weights=weights, ...)
      cv[i] <- sum(log(lamx)) - integral.im(lam)
    }
  }
  result <- bw.optim(cv, sigma, iopt=which.max(cv), 
                     creator="bw.ppl",
                     criterion="Likelihood Cross-Validation",
                     unitname=unitname(X))
  return(result)
}


bw.lppl <- function(X, ..., srange=NULL, ns=16, sigma=NULL,
                   weights=NULL, distance="euclidean", shortcut=FALSE) {
  stopifnot(is.lpp(X))
  if(!is.null(sigma)) {
    stopifnot(is.numeric(sigma) && is.vector(sigma))
    ns <- length(sigma)
  } else {
    if(!is.null(srange)) check.range(srange) else {
      dd <- diameter(Frame(X))                                               
      ss <- bw.scott.iso(X)
      srange <- range(c(ss/10, ss*5, dd/5))
    }
    sigma <- geomseq(from=srange[1L], to=srange[2L], length.out=ns)
  }
  cv <- numeric(ns)
  if(shortcut) {
    #' omit calculation of integral term
    #' precompute the geometry data
    lam1 <- density(X, sigma=sigma[1], weights=weights, distance=distance, ...,
                     savecomputed=TRUE)
    precooked <- attr(lam1, "savedstuff")
    for(i in 1:ns) {
      si <- sigma[i]
      lamx <- density(X, sigma=si, at="points", leaveoneout=TRUE,
                      weights=weights, distance=distance,
                      precomputed=precooked, ...)
      lamx <- pmax(0, lamx)
      cv[i] <- sum(log(lamx))
    }
  } else {
    #' full calculation
    precooked <- NULL
    cooking <- TRUE
    for(i in 1:ns) {
      si <- sigma[i]
      lamx <- density(X, sigma=si, at="points", leaveoneout=TRUE,
                      weights=weights, distance=distance, 
                      precomputed=precooked,
                      ...)
      lam <- density(X, sigma=si,
                     weights=weights, distance=distance,
                     precomputed=precooked,
                     savecomputed=cooking,
                     ...)
      if(cooking) {
        #' save geometry info for re-use in subsequent iterations
        precooked <- attr(lam, "savedstuff")
        cooking <- FALSE
      }
      lamx <- pmax(0, lamx)
      cv[i] <- sum(log(lamx)) - integral(lam)
    }
  }
  result <- bw.optim(cv, sigma, iopt=which.max(cv), 
                     creator="bw.lppl",
                     criterion="Likelihood Cross-Validation",
                     unitname=unitname(X))
  return(result)
}
