#'
#'     weightedStats.R
#'
#'   weighted versions of hist, var, median, quantile
#'
#'  $Revision: 1.2 $  $Date: 2016/07/17 05:26:54 $
#'


#'
#'    whist      weighted histogram
#'

whist <- function(x, breaks, weights=NULL) {
    N <- length(breaks)
    if(length(x) == 0) 
      h <- numeric(N+1)
    else {
      # classify data into histogram cells (breaks need not span range of data)
      cell <- findInterval(x, breaks, rightmost.closed=TRUE)
      # values of 'cell' range from 0 to N.
      nb <- N + 1L
      if(is.null(weights)) {
        ## histogram
        h <- tabulate(cell+1L, nbins=nb)
      } else {
        ##  weighted histogram
        if(!spatstat.options("Cwhist")) {
          cell <- factor(cell, levels=0:N)
          h <- unlist(lapply(split(weights, cell), sum, na.rm=TRUE))
        } else {
          h <- .Call("Cwhist",
                     as.integer(cell), as.double(weights), as.integer(nb),
                     PACKAGE = "spatstat")
        }
      }
    }
    h <- as.numeric(h)
    y <- h[2:N]
    attr(y, "low") <- h[1]
    attr(y, "high") <- h[N+1]
    return(y)
}

#' wrapper for computing weighted variance of a vector
#' Note: this includes a factor 1 - sum(v^2) in the denominator
#' where v = w/sum(w). See help(cov.wt)

weighted.var <- function(x, w, na.rm=TRUE) {
  bad <- is.na(w) | is.na(x)
  if(any(bad)) {
    if(!na.rm) return(NA_real_)
    ok <- !bad
    x <- x[ok]
    w <- w[ok]
  }
  cov.wt(matrix(x, ncol=1),w)$cov[]
}

#' weighted median

weighted.median <- function(x, w, na.rm=TRUE) {
  unname(weighted.quantile(x, probs=0.5, w=w, na.rm=na.rm))
}

#' weighted quantile

weighted.quantile <- function(x, w, probs=seq(0,1,0.25), na.rm=TRUE) {
  x <- as.numeric(as.vector(x))
  w <- as.numeric(as.vector(w))
  if(anyNA(x) || anyNA(w)) {
    ok <- !(is.na(x) | is.na(w))
    x <- x[ok]
    w <- w[ok]
  }
  stopifnot(all(w >= 0))
  if(all(w == 0)) stop("All weights are zero", call.=FALSE)
  #'
  oo <- order(x)
  x <- x[oo]
  w <- w[oo]
  Fx <- cumsum(w)/sum(w)
  #'
  result <- numeric(length(probs))
  for(i in seq_along(result)) {
    p <- probs[i]
    lefties <- which(Fx <= p)
    if(length(lefties) == 0) {
      result[i] <- x[1]
    } else {
      left <- max(lefties)
      result[i] <- x[left]
      if(Fx[left] < p && left < length(x)) {
        right <- left+1
        y <- x[left] + (x[right]-x[left]) * (p-Fx[left])/(Fx[right]-Fx[left])
        if(is.finite(y)) result[i] <- y
      }
    }
  }
  names(result) <- paste0(format(100 * probs, trim = TRUE), "%")
  return(result)
}

