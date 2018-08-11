#'
#'       matrixpower.R
#'
#'   $Revision: 1.1 $  $Date: 2016/11/13 01:50:51 $
#'

matrixsqrt <- function(x, complexOK=TRUE) {
  ## matrix square root
  if(length(dim(x)) != 2)
    stop("x must be a matrix")
  if(!is.matrix(x))
    x <- as.matrix(x)
  if(missing(complexOK) && is.complex(x)) complexOK <- TRUE
  if(!complexOK) stopifnot(is.numeric(x)) else
                 stopifnot(is.numeric(x) || is.complex(x))
  e <- eigen(x)
  values <- e$values
  vectors <- e$vectors
  if(any(values < 0)) {
    if(complexOK) values <- as.complex(values) else
    stop("matrix has negative eigenvalues: square root is complex",
         call.=FALSE)
  }
  y <- vectors %*% diag(sqrt(values)) %*% t(vectors)
  if(!is.null(dn <- dimnames(x)))
    dimnames(y) <- rev(dn)
  return(y)
}

matrixinvsqrt <- function(x, complexOK=TRUE) {
  ## matrix inverse square root
  if(length(dim(x)) != 2)
    stop("x must be a matrix")
  if(!is.matrix(x))
    x <- as.matrix(x)
  if(missing(complexOK) && is.complex(x)) complexOK <- TRUE
  if(!complexOK) stopifnot(is.numeric(x)) else
                 stopifnot(is.numeric(x) || is.complex(x))
  e <- eigen(x)
  values <- e$values
  vectors <- e$vectors
  if(any(values == 0))
    stop("matrix is singular; cannot compute inverse square root", call.=FALSE)
  if(any(values < 0)) {
    if(complexOK) values <- as.complex(values) else
    stop("matrix has negative eigenvalues: inverse square root is complex",
         call.=FALSE)
  }
  y <- vectors %*% diag(1/sqrt(values)) %*% t(vectors)
  if(!is.null(dn <- dimnames(x)))
    dimnames(y) <- rev(dn)
  return(y)
}

matrixpower <- function(x, power, complexOK=TRUE) {
  check.1.real(power)
  if(length(dim(x)) != 2)
    stop("x must be a matrix")
  if(!is.matrix(x))
    x <- as.matrix(x)
  if(missing(complexOK) && is.complex(x)) complexOK <- TRUE
  if(!complexOK) stopifnot(is.numeric(x)) else
                 stopifnot(is.numeric(x) || is.complex(x))
  e <- eigen(x)
  values <- e$values
  vectors <- e$vectors
  if(any(values == 0) && power < 0)
    stop("matrix is singular; cannot compute negative power", call.=FALSE)
  if(any(values < 0) && (power != ceiling(power))) {
    if(complexOK) values <- as.complex(values) else
    stop("matrix has negative eigenvalues: result is complex",
         call.=FALSE)
  }
  y <- vectors %*% diag(values^power) %*% t(vectors)
  if(!is.null(dn <- dimnames(x)))
    dimnames(y) <- rev(dn)
  return(y)
}
