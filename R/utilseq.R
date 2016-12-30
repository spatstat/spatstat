#'
#'   utilseq.R
#'
#'  Utilities for sequences, vectors, ranges of values
#'
#'       $Revision$ $Date$
#'

dropifsingle <- function(x) if(length(x) == 1) x[[1]] else x

#' ...............  ordering ......................

# defines the current favorite algorithm for 'order' 
fave.order <- function(x) { sort.list(x, method="quick", na.last=NA) }

# order statistic (for use in lapply calls) 
orderstats <- function(x, k, decreasing=FALSE) {
  if(decreasing) sort(x, decreasing=TRUE, na.last=TRUE)[k] else sort(x)[k]
}

# which value is k-th smallest
orderwhich <- function(x, k, decreasing=FALSE) {
  if(decreasing) order(x, decreasing=TRUE, na.last=TRUE)[k] else order(x)[k]
}


## ................ reverse cumulative sum .....................

revcumsum <- function(x) {
  #' equivalent to rev(cumsum(rev(x)))
  n <- length(x)
  if(identical(storage.mode(x), "integer")) {
    z <- .C("irevcumsum",
            x=as.integer(x),
            as.integer(n))
    return(z$x)
  } else {
    z <- .C("drevcumsum",
            x=as.double(x),
            as.integer(n))
    return(z$x)
  }
}

## ............. vectors of length 2 .................

as2vector <- function(x) {
  ## convert various wacky formats to numeric vector of length 2
  ## for use as coordinates of a single point.
  xname <- deparse(substitute(x))
  if(is.numeric(x)) {
    if(length(x) != 2)
      stop(paste(xname, "should have length 2"))
    return(x)
  }
  if(inherits(x, "ppp")) {
    #' coded so that it works even if spatstat is not loaded
    if(x$n != 1)
      stop(paste(xname, "should consist of exactly one point"))
    return(c(x$x, x$y))
  }
  if(is.list(x) && all(c("x", "y") %in% names(x))) {
    if(length(x$x) != 1) stop(paste0(xname, "$x should have length 1"))
    if(length(x$y) != 1) stop(paste0(xname, "$y should have length 1"))
    return(c(x$x, x$y))
  }
  stop(paste("Format of", sQuote(xname), "not understood"))
}

ensure2vector <- function(x) {
  xname <- deparse(substitute(x))
  if(!is.numeric(x))
    stop(paste(xname, "is not numeric"))
  n <- length(x)
  if(n == 0 || n > 2)
    stop(paste(xname, "should be of length 1 or 2"))
  if(n == 1)
    return(rep(x,2))
  return(x)
}


## ............. sequences ...................

prolongseq <- function(x, newrange, step=NULL) {
  ## Extend a sequence x so that it covers the new range.
  stopifnot(length(newrange) == 2 && newrange[1] < newrange[2])
  ## Check 'x' is an evenly-spaced sequence
  if(length(x) > 1) {
    dx <- diff(x)
    if(any(dx <= 0))
      stop("x must be an increasing sequence")
    if(diff(range(dx)) > 0.01 * abs(mean(dx)))
      stop("x must be evenly spaced")
  }
  ## Infer step length
  if(!is.null(step)) {
    check.1.real(step)
    stopifnot(step > 0)
  } else if(length(x) > 1) {
    step <- mean(dx)
  } else stop("step is needed when x is a single value")

  ## 
  if(max(x) < newrange[1] || min(x) > newrange[2])
    stop("x lies entirely outside the desired range")
    
  ## add or trim data to left
  if(x[1] > newrange[1]) {
    leftbit <- seq(from=x[1], to=newrange[1], by= -step)[-1]
    x <- c(rev(leftbit), x)
    nleft <- length(leftbit)
  } else {
    nx <- length(x)
    x <- x[x >= newrange[1]]
    nleft <- length(x) - nx
  }

  # add or trim data to right
  nx <- length(x)
  if(newrange[2] > x[nx]) {
    rightbit <- seq(from=x[nx], to=newrange[2], by= step)[-1]
    x <- c(x, rightbit)
    nright <- length(rightbit)
  } else {
    x <- x[x <= newrange[2]]
    nright <- length(x) - nx
  }
  attr(x, "nleft") <- nleft
  attr(x, "nright") <- nright
  return(x)
}

## fill gaps in a sequence
fillseq <- function(x, step=NULL) {
  xname <- short.deparse(substitute(x))
  n <- length(x)
  if(n <= 1) return(x)
  rx <- range(x)
  dx <- diff(x)
  if(any(dx < 0)) stop(paste(xname, "should be an increasing sequence"),
                       call.=FALSE)
  ## guess step length
  if(is.null(step)) {
    eps <- diff(rx)/1e7
    step <- min(dx[dx > eps])
  }
  ## make new sequence
  y <- seq(rx[1], rx[2], by=step)
  ny <- length(y)
  ## mapping from x to y
  i <- round((x - rx[1])/step) + 1L
  i <- pmin(ny, pmax(1, i))
  return(list(xnew=y, i=i))
}

geomseq <- function(from, to, length.out) {
  if(from <= 0 || to <= 0) stop("range limits must be positive")
  y <- exp(seq(from=log(from), to=log(to), length.out=length.out))
  y[1] <- from  #' avoid numerical error
  y[length.out] <- to
  return(y)
}

## ............. ranges ...................

intersect.ranges <- function(a, b, fatal=TRUE) {
  if(!is.null(a) && !is.null(b)) {
    lo <- max(a[1],b[1])
    hi <- min(a[2],b[2])
    if(lo <= hi)
      return(c(lo, hi))
  }
  if(fatal) stop("Intersection is empty")
  return(NULL)
}

inside.range <- function(x, r) {
  stopifnot(length(r) == 2 && r[1] <= r[2])
  return(x >= r[1] & x <= r[2])
}

check.in.range <- function(x, r, fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(inside.range(x, r))
    return(TRUE)
  if(fatal) 
    stop(paste(xname, "should be a number between",
               r[1], "and", r[2]),
         call.=FALSE)
  return(FALSE)
}

startinrange <- function(x0, dx, r) {
  ## find y = x0 + n * dx such that y \in r
  if(all(inside.range(x0, r))) return(x0)
  stopifnot(is.numeric(dx) && length(dx) == 1)
  y <- x0 + dx * round((mean(r) - x0)/dx)
  y[!inside.range(y, r)] <- NA
  return(y)
}

prettyinside <- function(x, ...) {
  r <- range(x, na.rm=TRUE)
  if(diff(r) == 0) return(r[1])
  p <- pretty(x, ...)
  ok <- inside.range(p, r)
  return(p[ok])
}

prettydiscrete <- function(x, n=10) {
  nx <- length(x)
  dx <- nx %/% n
  if(dx < 1) return(x)
  i <- 1 + (0:(n-1)) * dx
  return(x[i])
}


check.range <- function(x, fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(is.numeric(x) && identical(x, range(x, na.rm=TRUE)))
    return(TRUE)
  if(fatal) 
    stop(paste(xname, "should be a vector of length 2 giving (min, max)"))
  return(FALSE)
}

evenly.spaced <- function(x, tol=1e-07) {
  # test whether x is evenly spaced and increasing
  dx <- diff(x)
  if(any(dx <= .Machine$double.eps))
    return(FALSE)
  # The following test for equal spacing is used in hist.default
  if(diff(range(dx)) > tol * mean(dx))
    return(FALSE)
  return(TRUE)
}

equispaced <- function(z, reltol=0.001) {
  dz <- diff(z)
  return(diff(range(dz)) < reltol * mean(dz))
}


adjustthinrange <- function(ur,vstep,vr) {
  if(diff(ur) >= vstep) return(ur)
  ur <- mean(ur) + c(-1,1) * vstep/2
  if(ur[1] < vr[1]) ur <- vr[1] + c(0,1)*vstep
  if(ur[2] > vr[2]) ur <- vr[2] - c(1,0)*vstep
  return(ur)
}

fastFindInterval <- function(x, b, labels=FALSE, reltol=0.001) {
  nintervals <- length(b) - 1
  nx <- length(x)
  if(nx == 0)
    return(rep(0, nintervals))
  ##
  if(equispaced(b, reltol)) {
    ## breaks are equally spaced
    zz <- .C("fastinterv",
             x          = as.double(x),
             n          = as.integer(nx),
             brange     = as.double(range(b)),
             nintervals = as.integer(nintervals),
             y          = as.integer(integer(nx))
             )
    y <- zz$y
  } else {
    ## use R's interval search algorithm
    y <- findInterval(x, b, rightmost.closed=TRUE)
  }
  if(labels) {
    blab <- paste0("[",
                   b[1:nintervals],
                   ",",
                   b[-1],
                   c(rep(")", nintervals-1), "]"))
    y <- as.integer(y)
    levels(y) <- as.character(blab)
    class(y) <- "factor"
  }
  return(y)
}

# ...................................................
# efficient replacements for ifelse()
# 'a' and 'b' are single values
# 'x' and 'y' are vectors of the same length as 'test'

# ifelse(test, a, b)
ifelseAB <- function(test,  a, b) {
  y <- rep.int(b, length(test))
  y[test] <- a
  return(y)
}

# ifelse(test, a, x)
ifelseAX <- function(test, a, x) {
  y <- x
  y[test] <- a
  return(y)
}

# ifelse(test, x, b)
ifelseXB <- function(test, x, b) {
  y <- rep.int(b, length(test))
  y[test] <- x[test]
  return(y)
}
  
# ifelse(test, x, y)
ifelseXY <- function(test, x, y) {
  z <- y
  z[test] <- x[test]
  return(z)
}

#.... very special cases ......

# ifelse(test, 1, NA)
ifelse1NA <- function(test) {
  y <- as.integer(test)
  y[!test] <- NA
  return(y)
}

# ifelse(test, 0, NA)
ifelse0NA <- function(test) {
  nyet <- !test
  y <- as.integer(nyet)
  y[nyet] <- NA
  return(y)
}

# ifelse(test, -x, x)
ifelseNegPos <- function(test, x) {
  y <- x
  y[test] <- -x[test]
  return(y)
}


ratiotweak <- function(a, b, overzero=NA, zerozero=overzero) {
  # map x/0 to 'overzero' and 0/0 to 'zerozero'
  result <- a/b
  bzero <- (b == 0)
  result[ bzero ] <- overzero
  if(!missing(zerozero)) {
    abzero <- bzero & (a == 0)
    result[ abzero ] <- zerozero
  }
  return(result)
}

natozero <- function(x) {
  #' map NA to zero (e.g. in tapply)
  x[is.na(x)] <- 0
  return(x)
}

insertinlist <- function(x, i, y) {
  ## insert a possibly longer or shorter list 'y'
  ## into serial position 'i' in list 'x'
  n <- length(x)
  if(n == 0) return(y)
  m <- seq_len(n)
  names(m) <- names(x)
  i <- m[[i]] # convert 'i' to integer index
  stopifnot(length(i) == 1)
  if(n == 1) return(y)
  xleft <- x[seq_len(i-1)]
  xright <- x[i + seq_len(n-i)]
  z <- c(xleft, y, xright)
  return(z)
}

#' ..... rounding ..............................

dround <- function(x) {
  round(x, getOption('digits'))
}

niceround <- function(x, m=c(1,2,5,10)) {
  expo <- 10^as.integer(floor(log10(x)))
  y <- m * expo
  z <- y[which.min(abs(y - x))]
  return(z)
}


