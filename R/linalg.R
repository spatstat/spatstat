#
# linalg.R
#
# $Revision: 1.4 $ $Date: 2011/08/05 01:47:50 $
#

sumouter <- function(x, w) {
  stopifnot(is.matrix(x))
  p <- ncol(x)
  n <- nrow(x)
  nama <- colnames(x)
  # transpose (compute outer squares of columns)
  tx <- t(x)
  ok <- apply(is.finite(tx), 2, all)
  if(missing(w)) {
    w <- rep(1, n)
  } else {
    if(length(w) != n)
      stop(paste("The length of w does not match the number of rows of x",
                 "\t(", length(w), "!=", n, ")"))
    ok <- ok & is.finite(w)
  }
  if(!all(ok)) {
    tx <- tx[ , ok, drop=FALSE]
    w <- w[ok]
  }
  DUP <- spatstat.options("dupC")
  z <- .C("wsumouter",
          x=as.double(tx),
          n=as.integer(n),
          p=as.integer(p),
          w=as.double(w),
          y=as.double(numeric(p * p)),
          DUP=DUP,
          PACKAGE="spatstat")
  out <- matrix(z$y, p, p)
  if(!is.null(nama))
     dimnames(out) <- list(nama, nama)
  return(out)
}

quadform <- function(x, v) {
  stopifnot(is.matrix(x))
  p <- ncol(x)
  n <- nrow(x)
  nama <- rownames(x)
  # transpose (evaluate quadratic form for each column)
  tx <- t(x)
  ok <- apply(is.finite(tx), 2, all)
  allok <- all(ok)
  if(!allok) {
    tx <- tx[ , ok, drop=FALSE]
    n <- ncol(tx)
  }
  if(missing(v)) {
    v <- diag(rep(1, p))
  } else {
    stopifnot(is.matrix(v))
    if(nrow(v) != ncol(v)) stop("v should be a square matrix")
    stopifnot(ncol(x) == nrow(v))
  }
  DUP <- spatstat.options("dupC")  
  z <- .C("quadform",
          x=as.double(tx),
          n=as.integer(n),
          p=as.integer(p),
          v=as.double(v),
          y=as.double(numeric(n)),
          DUP=DUP,
          PACKAGE="spatstat")
  result <- z$y
  names(result) <- nama[ok]
  if(allok)
    return(result)
  fullresult <- rep(NA, length(ok))
  fullresult[ok] <- result
  names(fullresult) <- nama
  return(fullresult)
}


