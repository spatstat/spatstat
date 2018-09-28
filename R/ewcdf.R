#
#     ewcdf.R
#
#     $Revision: 1.13 $  $Date: 2018/09/28 05:10:02 $
#
#  With contributions from Kevin Ummel
#

ewcdf <- function(x, weights=rep(1/length(x), length(x)))
{
  nw <- length(weights)
  nx <- length(x)
  if(nw == 0) {
    weights <- rep(1/nx, nx)
  } else if(nw == 1) {
    weights <- rep(weights, nx)
  } else if(nw != nx) stopifnot(length(x) == length(weights))
  # remove NA's together
  nbg <- is.na(x) 
  x <- x[!nbg]
  weights <- weights[!nbg]
  n <- length(x)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")
  stopifnot(all(weights >= 0))
  # sort in increasing order of x value
  ox <- fave.order(x)
  x <- x[ox]
  w <- weights[ox]
  # find jump locations and match
  vals <- sortunique(x)
  xmatch <- factor(match(x, vals), levels=seq_along(vals))
  # sum weight in each interval
  wmatch <- tapply(w, xmatch, sum)
  wmatch[is.na(wmatch)] <- 0
  cumwt <- cumsum(wmatch)
  # make function
  rval <- approxfun(vals, cumwt,
                    method = "constant", yleft = 0, yright = sum(wmatch),
                    f = 0, ties = "ordered")
  class(rval) <- c("ewcdf", "ecdf", "stepfun", class(rval))
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

