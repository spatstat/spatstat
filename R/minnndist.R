##
##  minnndist.R
##
## Fast versions of min(nndist(X)), max(nndist(X))
##
##  $Revision: 1.5 $  $Date: 2017/06/05 10:31:58 $

minnndist <- function(X, positive=FALSE) {
  stopifnot(is.ppp(X))
  n <- npoints(X)
  if(n <= 1) return(NA)
  x <- X$x
  y <- X$y
  o <- fave.order(y)
  big <- sqrt(.Machine$double.xmax)
  if(positive) {
      z <- .C("minPnnd2",
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              PACKAGE = "spatstat")
  } else {
      z <- .C("minnnd2",
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              PACKAGE = "spatstat")
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
  big <- sqrt(.Machine$double.xmax)
  if(positive) {
      z <- .C("maxPnnd2",
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              PACKAGE = "spatstat")
  } else {
      z <- .C("maxnnd2",
              n = as.integer(n),
              x = as.double(x[o]),
              y = as.double(y[o]),
              as.double(big),
              result = as.double(numeric(1)),
              PACKAGE = "spatstat")
  }
  return(sqrt(z$result))
}

