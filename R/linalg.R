#
# linalg.R
#
# $Revision: 1.6 $ $Date: 2013/04/25 06:37:43 $
#

sumouter <- function(x, w=NULL) {
  stopifnot(is.matrix(x))
  p <- ncol(x)
  n <- nrow(x)
  nama <- colnames(x)
  # transpose (compute outer squares of columns)
  tx <- t(x)
  ok <- apply(is.finite(tx), 2, all)
  if(!is.null(w)) {
    if(length(w) != n)
      stop(paste("The length of w does not match the number of rows of x",
                 "\t(", length(w), "!=", n, ")"))
    ok <- ok & is.finite(w)
  }
  if(!all(ok)) {
    tx <- tx[ , ok, drop=FALSE]
    if(!is.null(w)) w <- w[ok]
  }
  DUP <- spatstat.options("dupC")
  if(is.null(w)) {
    z <- .C("sumouter",
            x=as.double(tx),
            n=as.integer(n),
            p=as.integer(p),
            y=as.double(numeric(p * p)),
            DUP=DUP,
            PACKAGE="spatstat")
  } else {
    z <- .C("wsumouter",
            x=as.double(tx),
            n=as.integer(n),
            p=as.integer(p),
            w=as.double(w),
            y=as.double(numeric(p * p)),
            DUP=DUP,
            PACKAGE="spatstat")
  }
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
    v <- diag(rep.int(1, p))
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
  fullresult <- rep.int(NA_real_, length(ok))
  fullresult[ok] <- result
  names(fullresult) <- nama
  return(fullresult)
}

sumsymouter <- function(x, w=NULL) {
  # computes the sum of outer(x[,i,j], x[,j,i]) * w[i,j] over all pairs i != j
  stopifnot(is.array(x) && length(dim(x)) == 3)
  if(dim(x)[2] != dim(x)[3])
    stop("The second and third dimensions of x should be equal")
  if(!is.null(w)) {
    stopifnot(is.matrix(w))
    if(!all(dim(w) == dim(x)[-1]))
      stop("Dimensions of w should match the second and third dimensions of x")
  }
  p <- dim(x)[1]
  n <- dim(x)[2]
  DUP <- spatstat.options("dupC")
  if(is.null(w)) {
    zz <- .C("sumsymouter",
             x = as.double(x),
             p = as.integer(p),
             n = as.integer(n),
             y = as.double(numeric(p * p)),
             DUP = DUP,
             PACKAGE = "spatstat")
  } else {
    zz <- .C("wsumsymouter",
             x = as.double(x),
             w = as.double(w),
             p = as.integer(p),
             n = as.integer(n),
             y = as.double(numeric(p * p)),
             DUP = DUP,
             PACKAGE = "spatstat")
  }
  matrix(zz$y, p, p)
}

