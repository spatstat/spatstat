#
#     ewcdf.R
#
#     $Revision: 1.18 $  $Date: 2019/03/12 11:14:36 $
#
#  With contributions from Kevin Ummel
#

ewcdf <- function(x, weights=NULL, normalise=TRUE, adjust=1)
{
  nx <- length(x)
  nw <- length(weights)
  weighted <- (nw > 0)

  if(weighted) {
    check.nvector(weights, things="entries of x", oneok=TRUE)
    stopifnot(all(weights >= 0))
    if(nw == 1) 
      weights <- rep(weights, nx)
  }
  
  ## remove NA's
  nbg <- is.na(x) 
  x <- x[!nbg]
  if(weighted) weights <- weights[!nbg]
  n <- length(x)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")

  ## sort in increasing order of x value
  if(!weighted) {
    x <- sort(x)
    w <- rep(1, n)
  } else {
    ox <- fave.order(x)
    x <- x[ox]
    w <- weights[ox]
  }
  ## find jump locations and match
  rl <- rle(x)
  vals <- rl$values
  if(!weighted) {
    wmatch <- rl$lengths
  } else {
    nv <- length(vals)
    wmatch <- .C("tabsumweight",
                 nx=as.integer(n),
                 x=as.double(x),
                 w=as.double(w),
                 nv=as.integer(nv),
                 v=as.double(vals),
                 z=as.double(numeric(nv)),
                 PACKAGE="spatstat")$z
  }
  ## cumulative weight in each interval
  cumwt <- cumsum(wmatch)
  totwt <- sum(wmatch)
  ## rescale ?
  if(normalise) {
    cumwt <- cumwt/totwt
    totwt <- 1
  } else if(adjust != 1) {
    cumwt <- adjust * cumwt
    totwt <- adjust * totwt
  }
  ## make function
  rval <- approxfun(vals, cumwt,
                    method = "constant", yleft = 0, yright = totwt,
                    f = 0, ties = "ordered")
  class(rval) <- c("ewcdf",
                   if(normalise) "ecdf" else NULL,
                   "stepfun", class(rval))
  assign("w", w, envir=environment(rval))
  attr(rval, "call") <- sys.call()
  return(rval)
}

  # Hacked from stats:::print.ecdf
print.ewcdf <- function (x, digits = getOption("digits") - 2L, ...) {
  cat("Weighted empirical CDF \nCall: ")
  print(attr(x, "call"), ...)
  env <- environment(x)
  xx <- get("x", envir=env)
  ww <- get("w", envir=env)
  n <- length(xx)
  i1 <- 1L:min(3L, n)
  i2 <- if (n >= 4L) max(4L, n - 1L):n else integer()
  numform <- function(x) paste(formatC(x, digits = digits), collapse = ", ")
  cat(" x[1:", n, "] = ", numform(xx[i1]), if (n > 3L) 
      ", ", if (n > 5L) 
      " ..., ", numform(xx[i2]), "\n", sep = "")
  cat(" weights[1:", n, "] = ", numform(ww[i1]), if (n > 3L) 
      ", ", if (n > 5L) 
      " ..., ", numform(ww[i2]), "\n", sep = "")
  invisible(x)
}

quantile.ewcdf <- function(x, probs=seq(0,1,0.25), names=TRUE, ...,
                           normalise=TRUE, type=1) {
  trap.extra.arguments(..., .Context="quantile.ewcdf")
  if(!(type %in% c(1,2)))
    stop("Only quantiles of type 1 and 2 are implemented", call.=FALSE)
  env <- environment(x)
  xx <- get("x", envir=env)
  n <- length(xx)
  Fxx <- get("y", envir=env)
  maxFxx <- max(Fxx)
  eps <- 100 * .Machine$double.eps
  if(normalise) {
    Fxx <- Fxx/maxFxx
    maxp <- 1
  } else {
    maxp <- maxFxx
  }
  if(any((p.ok <- !is.na(probs)) &
         (probs/maxp < -eps | probs/maxp > 1 + eps))) {
    allowed <- if(normalise) "[0,1]" else
               paste("permitted range", prange(c(0, maxp)))
    stop(paste("'probs' outside", allowed), call.=FALSE)
  }
  if (na.p <- any(!p.ok)) {
    o.pr <- probs
    probs <- probs[p.ok]
    probs <- pmax(0, pmin(maxp, probs))
  }
  np <- length(probs)
  if (n > 0 && np > 0) {
    qs <- numeric(np)
    if(type == 1) {
      ## right-continuous inverse
      for(k in 1:np) qs[k] <- xx[min(which(Fxx >= probs[k]))]
    } else {
      ## average of left and right continuous
      for(k in 1:np) {
        pk <- probs[k]
        ik <- min(which(Fxx >= probs[k]))
        qs[k] <- if(Fxx[ik] > pk) (xx[ik] + xx[ik-1L])/2 else xx[ik]
      }
    }
  } else {
    qs <- rep(NA_real_, np)
  }
  if (names && np > 0L) {
    dig <- max(2L, getOption("digits"))
    if(normalise) {
      probnames <-
        if(np < 100) formatC(100 * probs, format="fg", width=1, digits=dig) else
        format(100 * probs, trim = TRUE, digits = dig)
      names(qs) <- paste0(probnames, "%")
    } else {
      names(qs) <-
        if(np < 100) formatC(probs, format="fg", width=1, digits=dig) else
        format(probs, trim=TRUE, digits=dig)
    }
  }
  if (na.p) {
    o.pr[p.ok] <- qs
    names(o.pr) <- rep("", length(o.pr))
    names(o.pr)[p.ok] <- names(qs)
    o.pr
  } else qs
}

