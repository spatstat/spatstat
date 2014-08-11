##
##  minnndist.R
##
## Fast versions of min(nndist(X)), max(nndist(X))
##
##  $Revision: 1.2 $  $Date: 2014/03/24 07:31:07 $

minnndist <- function(X) {
  stopifnot(is.ppp(X))
  n <- npoints(X)
  if(n <= 1) return(NA)
  x <- X$x
  y <- X$y
  o <- fave.order(y)
  DUP <- spatstat.options("dupC")
  big <- sqrt(.Machine$double.xmax)
  z <- .C("minnnd2",
          n = as.integer(n),
          x = as.double(x[o]),
          y = as.double(y[o]),
          as.double(big),
          result = as.double(numeric(1)),
          DUP=DUP)
  return(sqrt(z$result))
}

maxnndist <- function(X) {
  stopifnot(is.ppp(X))
  n <- npoints(X)
  if(n <= 1) return(NA)
  x <- X$x
  y <- X$y
  o <- fave.order(y)
  DUP <- spatstat.options("dupC")
  big <- sqrt(.Machine$double.xmax)
  z <- .C("maxnnd2",
          n = as.integer(n),
          x = as.double(x[o]),
          y = as.double(y[o]),
          as.double(big),
          result = as.double(numeric(1)),
          DUP=DUP)
  return(sqrt(z$result))
}

          
