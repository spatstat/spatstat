##
##  minnndist.R
##
## Fast versions of min(nndist(X)), max(nndist(X))
##
##  $Revision: 1.3 $  $Date: 2014/09/18 01:17:42 $

minnndist <- function(X, positive=FALSE) {
  stopifnot(is.ppp(X))
  n <- npoints(X)
  if(n <= 1) return(NA)
  x <- X$x
  y <- X$y
  o <- fave.order(y)
  DUP <- spatstat.options("dupC")
  big <- sqrt(.Machine$double.xmax)
  if(positive) {
      z <- .C("minPnnd2",
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              DUP=DUP)
  } else {
      z <- .C("minnnd2",
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              DUP=DUP)
  }
  return(sqrt(z$result))
}

maxnndist <- function(X, positive=FALSE) {
  stopifnot(is.ppp(X))
  n <- npoints(X)
  if(n <= 1) return(NA)
  x <- X$x
  y <- X$y
  o <- fave.order(y)
  DUP <- spatstat.options("dupC")
  big <- sqrt(.Machine$double.xmax)
  if(positive) {
      z <- .C("maxPnnd2",
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              DUP=DUP)
  } else {
      z <- .C("maxnnd2",
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              DUP=DUP)
  }
  return(sqrt(z$result))
}

          
