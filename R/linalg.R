#
# linalg.R
#
#  Linear Algebra
#
# $Revision: 1.31 $ $Date: 2020/05/10 00:59:17 $
#

sumouter <- function(x, w=NULL, y=x) {
  #' compute matrix sum_i (w[i] * outer(x[i,], y[i,]))
  stopifnot(is.matrix(x))
  weighted <- !is.null(w)
  symmetric <- missing(y) || identical(x,y)
  if(!symmetric) {
    stopifnot(is.matrix(y))
    stopifnot(nrow(x) == nrow(y))
  }
  if(weighted) {
    if(max(1L, dim(w)) > 1L) stop("w should be a vector")
    w <- as.vector(w)
    if(is.complex(w)) {
      stopifnot(length(w) == nrow(x))
    } else {
      w <- as.numeric(w)
      check.nvector(w, nrow(x), things="rows of x")
    }
  }
    
  #' ............ complex arguments ....................
  #' handle complex weights ...
  if(weighted && is.complex(w)) {
    a <- if(symmetric) sumouter(x, Re(w)) else sumouter(x, Re(w), y)
    b <- if(symmetric) sumouter(x, Im(w)) else sumouter(x, Im(w), y)
    result <- a + b * 1i
    return(result)
  }
  #' handle complex x, y ...
  if(is.complex(x) || is.complex(y)) {
    Rex <- Re(x)
    Imx <- Im(x)
    Sux <- Rex + Imx
    if(symmetric) {
      RexRey <- sumouter(Rex, w)
      ImxImy <- sumouter(Imx, w)
      SuxSuy <- sumouter(Sux, w)
    } else {
      Rey <- Re(y)
      Imy <- Im(y)
      Suy <- Rey + Imy
      RexRey <- sumouter(Rex, w, Rey)
      ImxImy <- sumouter(Imx, w, Imy)
      SuxSuy <- sumouter(Sux, w, Suy)
    }
    result <- RexRey - ImxImy + (SuxSuy - RexRey - ImxImy) * 1i
    return(result)
  }
  #' .......... All arguments are numeric. .................
  #' Transpose (compute outer squares of columns)
  tx <- t(x)
  if(!symmetric) ty <- t(y)
  #' check for NA etc
  ok <- apply(is.finite(tx), 2, all)
  if(!symmetric) ok <- ok & apply(is.finite(ty), 2, all)
  if(weighted) ok <- ok & is.finite(w)
  #' remove NA etc
  if(!all(ok)) {
    tx <- tx[ , ok, drop=FALSE]
    if(!symmetric) ty <- ty[ , ok, drop=FALSE]
    if(weighted) w <- w[ok]
  }
  #' call C code
  if(symmetric) {
    n <- ncol(tx)
    p <- nrow(tx)
    if(is.null(w)) {
      zz <- .C("Csumouter",
               x=as.double(tx),
               n=as.integer(n),
               p=as.integer(p),
               y=as.double(numeric(p * p)),
               PACKAGE = "spatstat")
    } else {
      zz <- .C("Cwsumouter",
               x=as.double(tx),
               n=as.integer(n),
               p=as.integer(p),
               w=as.double(w),
               y=as.double(numeric(p * p)),
               PACKAGE = "spatstat")
    }
    out <- matrix(zz$y, p, p)
    if(!is.null(nama <- colnames(x)))
      dimnames(out) <- list(nama, nama)
  } else {
    n <- ncol(tx)
    px <- nrow(tx)
    py <- nrow(ty)
    if(is.null(w)) {
      zz <- .C("Csum2outer",
               x=as.double(tx),
               y=as.double(ty),
               n=as.integer(n),
               px=as.integer(px),
               py=as.integer(py),
               z=as.double(numeric(px * py)),
               PACKAGE = "spatstat")
    } else {
      zz <- .C("Cwsum2outer",
               x=as.double(tx),
               y=as.double(ty),
               n=as.integer(n),
               px=as.integer(px),
               py=as.integer(py),
               w=as.double(w),
               z=as.double(numeric(px * py)))
    }
    out <- matrix(zz$z, px, py)
    namx <- colnames(x)
    namy <- colnames(y)
    if(!is.null(namx) || !is.null(namy))
      dimnames(out) <- list(namx, namy)
  }
  return(out)
}

quadform <- function(x, v) {
  #' compute vector of values y[i] = x[i, ] %*% v %*% t(x[i,])
  stopifnot(is.matrix(x))
  #'
  if(missing(v)) {
    v <- diag(rep.int(1, p))
  } else {
    stopifnot(is.matrix(v))
    if(nrow(v) != ncol(v)) stop("v should be a square matrix")
    stopifnot(ncol(x) == nrow(v))
  }
  #' handle complex values
  if(is.complex(v)) {
    a <- quadform(x, Re(v))
    b <- quadform(x, Im(v))
    result <- a + b * 1i
    return(result)
  }
  if(is.complex(x)) {
    A <- quadform(Re(x), v)
    B <- quadform(Im(x), v)
    D <- quadform(Re(x) + Im(x), v)
    result <- A - B + (D - A - B) * 1i
    return(result)
  }
  #' arguments are numeric
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
  z <- .C("Cquadform",
          x=as.double(tx),
          n=as.integer(n),
          p=as.integer(p),
          v=as.double(v),
          y=as.double(numeric(n)),
          PACKAGE = "spatstat")
  result <- z$y
  names(result) <- nama[ok]
  if(allok)
    return(result)
  fullresult <- rep.int(NA_real_, length(ok))
  fullresult[ok] <- result
  names(fullresult) <- nama
  return(fullresult)
}

bilinearform <- function(x, v, y) {
  #' compute vector of values z[i] = x[i, ] %*% v %*% t(y[i,])
  stopifnot(is.matrix(x))
  stopifnot(is.matrix(y))
  stopifnot(identical(dim(x), dim(y)))
  #' handle complex values
  if(is.complex(v)) {
    a <- bilinearform(x, Re(v), y)
    b <- bilinearform(x, Im(v), y)
    result <- a + b * 1i
    return(result)
  }
  if(is.complex(x) || is.complex(y)) {
    a <- bilinearform(Re(x), v, Re(y))
    b <- bilinearform(Im(x), v, Im(y))
    d <- bilinearform(Re(x)+Im(x), v, Re(y)+Im(y))
    result <- a - b + (d - a - b) * 1i
    return(result)
  }
  #' arguments are numeric
  p <- ncol(x)
  n <- nrow(x)
  nama <- rownames(x)
  # transpose (evaluate quadratic form for each column)
  tx <- t(x)
  ty <- t(y)
  ok <- matcolall(is.finite(tx)) & matcolall(is.finite(ty))
  allok <- all(ok)
  if(!allok) {
    tx <- tx[ , ok, drop=FALSE]
    ty <- ty[ , ok, drop=FALSE]
    n <- ncol(tx)
  }
  if(missing(v)) {
    v <- diag(rep.int(1, p))
  } else {
    stopifnot(is.matrix(v))
    if(nrow(v) != ncol(v)) stop("v should be a square matrix")
    stopifnot(ncol(x) == nrow(v))
  }
  z <- .C("Cbiform",
          x=as.double(tx),
          y=as.double(ty),
          n=as.integer(n),
          p=as.integer(p),
          v=as.double(v),
          z=as.double(numeric(n)),
          PACKAGE = "spatstat")
  result <- z$z
  names(result) <- nama[ok]
  if(allok)
    return(result)
  fullresult <- rep.int(NA_real_, length(ok))
  fullresult[ok] <- result
  names(fullresult) <- nama
  return(fullresult)
}

sumsymouter <- function(x, w=NULL, distinct=TRUE) {
  ## x is a 3D array
  ## w is a matrix
  ## Computes the sum of outer(x[,i,j], x[,j,i]) * w[i,j] over all pairs i != j
  ## handle complex values
  if(is.complex(w)) {
    a <- sumsymouter(x, Re(w), distinct=distinct)
    b <- sumsymouter(x, Im(w), distinct=distinct)
    result <- a + b * 1i
    return(result)
  }
  if(is.complex(x)) {
    a <- sumsymouter(Re(x), w=w, distinct=distinct)
    b <- sumsymouter(Im(x), w=w, distinct=distinct)
    d <- sumsymouter(Re(x)+Im(x), w=w, distinct=distinct)
    result <- a - b + (d - a - b) * 1i
    return(result)
  }
  ## handle sparse arrays
  if(inherits(x, c("sparseSlab", "sparse3Darray")) &&
     (is.null(w) || inherits(w, "sparseMatrix")))
    return(sumsymouterSparse(x, w, distinct=distinct))
  ## arguments are numeric
  x <- as.array(x)
  stopifnot(length(dim(x)) == 3)
  if(dim(x)[2L] != dim(x)[3L])
    stop("The second and third dimensions of x should be equal")
  if(!is.null(w)) {
    w <- as.matrix(w)
    if(!all(dim(w) == dim(x)[-1L]))
      stop("Dimensions of w should match the second and third dimensions of x")
  }
  p <- dim(x)[1L]
  n <- dim(x)[2L]
  if(!distinct) {
    ## contributions from all pairs i,j
    if(is.null(w)) {
      zz <- .C("Csumsymouter",
               x = as.double(x),
               p = as.integer(p),
               n = as.integer(n),
               y = as.double(numeric(p * p)),
               PACKAGE = "spatstat")
    } else {
      zz <- .C("Cwsumsymouter",
               x = as.double(x),
               w = as.double(w),
               p = as.integer(p),
               n = as.integer(n),
               y = as.double(numeric(p * p)),
               PACKAGE = "spatstat")
    }
  } else {
    ## contributions from pairs i != j
    if(is.null(w)) {
      zz <- .C("CsumDsymouter",
               x = as.double(x),
               p = as.integer(p),
               n = as.integer(n),
               y = as.double(numeric(p * p)),
               PACKAGE = "spatstat")
    } else {
      zz <- .C("CwsumDsymouter",
               x = as.double(x),
               w = as.double(w),
               p = as.integer(p),
               n = as.integer(n),
               y = as.double(numeric(p * p)),
               PACKAGE = "spatstat")
    }
  }
  matrix(zz$y, p, p)
}

## matrix utilities

checksolve <- function(M, action=c("fatal", "warn", "silent"),
                       descrip, target="inverse") {
  Mname <- short.deparse(substitute(M))
  action <- match.arg(action)
  Minv <- try(solve(M), silent=(action=="silent"))
  if(!inherits(Minv, "try-error"))
    return(Minv)
  if(missing(descrip) || sum(nzchar(descrip)) == 0)
    descrip <- paste("the matrix", sQuote(Mname))
  whinge <- paste("Cannot compute", paste0(target, ":"),
                  descrip, "is singular")
  switch(action,
         fatal=stop(whinge, call.=FALSE),
         warn= warning(whinge, call.=FALSE),
         silent={})
  return(NULL)
}

check.mat.mul <- function(A, B, Acols="columns of A", Brows="rows of B",
                          fatal=TRUE) {
  # check whether A %*% B would be valid: if not, print a useful message
  if(!is.matrix(A)) A <- matrix(A, nrow=1, dimnames=list(NULL, names(A)))
  if(!is.matrix(B)) B <- matrix(B, ncol=1, dimnames=list(names(B), NULL))
  nA <- ncol(A)
  nB <- nrow(B) 
  if(nA == nB) return(TRUE)
  if(!fatal) return(FALSE)
  if(any(nzchar(Anames <- colnames(A))))
    message(paste0("Names of ", Acols, ": ", commasep(Anames)))
  if(any(nzchar(Bnames <- rownames(B))))
    message(paste0("Names of ", Brows, ": ", commasep(Bnames)))
  stop(paste("Internal error: number of", Acols, paren(nA),
             "does not match number of", Brows, paren(nB)),
       call.=FALSE)
}

