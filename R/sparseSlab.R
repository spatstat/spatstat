#'
#' sparseSlab.R
#'
#' Sparse 3D arrays or 'slabs' (where one of the dimensions is small)
#' represented as lists of sparse 2D matrices
#' 
#' $Revision: 1.5 $  $Date: 2016/02/18 07:01:21 $
#'

sparseSlab <- function(matlist, stackdim=3) {
  stopifnot(is.list(matlist))
  if(is.null(stackdim)) stackdim <- 3
  stopifnot(stackdim %in% 1:3)
  dims <- lapply(matlist, dim)
  if(any(lengths(dims) != 2))
    stop("Elements of matlist must be matrices [sparse or full]")
  if(length(unique(dims)) != 1)
    stop("Matrices must have identical dimensions")
  if(length(unique(lapply(matlist, class))) != 1)
    stop("Matrices must belong to the same class")
  d <- integer(3)
  d[stackdim] <- length(matlist)
  d[-stackdim] <- dims[[1]]
  result <- list(m=matlist, dim=d, stackdim=stackdim)
  class(result) <- "sparseSlab"
  return(result)
}

as.sparseSlab <- function(x, ..., stackdim=NULL) {
  if(inherits(x, "sparseSlab")) {
    y <- x
  } else if(inherits(x, "sparseMatrix")) {
    y <- sparseSlab(list(x), stackdim)
  } else if(is.matrix(x)) {
    y <- sparseSlab(list(as(x, "sparseMatrix")), stackdim)
  } else if(is.list(x) && all(sapply(x, is.matrix))) {
    y <- sparseSlab(lapply(x, as, Class="sparseMatrix"), stackdim)
  } else if(is.array(x)) {
    if(is.null(stackdim)) stackdim <- which.min(dim(x))
    stopifnot(stackdim %in% 1:3)
    n <- dim(x)[stackdim]
    d <- dim(x)[-stackdim]
    matlist <- vector(mode="list", length=n)
    for(i in seq_len(n)) {
      mi <- switch(stackdim,
                   matrix(x[i,,], d[1], d[2]),
                   matrix(x[,i,], d[1], d[2]),
                   matrix(x[,,i], d[1], d[2]))
      matlist[[i]] <- as(mi, Class="sparseMatrix")
    }
    y <- sparseSlab(matlist, stackdim)
  } else {          
    warning("I don't know how to convert x to a sparse array")
    return(NULL)
  }
  if(!is.null(stackdim) && stackdim != y$stackdim) {
    d <- integer(3)
    d[stackdim] <- dim(y)[y$stackdim]
    d[-stackdim] <- dim(y)[-y$stackdim]
    y$dim <- d
    y$stackdim <- stackdim
  }
  return(y)
}

dim.sparseSlab <- function(x) { x$dim }

dimnames.sparseSlab <- function(x) {
  with(x, {
    out <- vector(mode="list", length=3)
    out[stackdim] <- list(names(m))
    out[-stackdim] <- dimnames(m[[1]]) %orifnull% list(NULL, NULL)
    return(out)
  })
}

"dimnames<-.sparseSlab" <- function(x, value) {
  if(!is.list(value)) value <- list(value)
  if(length(value) == 1) value <- rep(value, 3)
  k <- x$stackdim
  names(x$m) <- value[[k]]
  othernames <- value[-k]
  for(i in seq_along(x$m))
    dimnames(x$m[[i]]) <- othernames
  return(x)
}

print.sparseSlab <- function(x, ...) {
  cat("Sparse array of dimensions", paste(dim(x), collapse="x"), fill=TRUE)
  stackdim <- x$stackdim
  m <- x$m
  cat("(sparse matrices stacked along dimension", paste0(stackdim, ")"),
      fill=TRUE)
  ns <- dim(x)[stackdim]
  indices <- character(3)
  z <- if(!is.null(names(m))) shQuote(names(m)) else as.character(seq_len(ns))
  for(k in seq_len(ns)) {
    indices[stackdim] <- z[k]
    cat(paste0("\n\t[", paste(indices, collapse=", "), "]\n\n"))
    print(m[[k]])
  }
  return(invisible(NULL))
}

aperm.sparseSlab <- function(a, perm=NULL, resize=TRUE, ...) {
  if(!resize)
    stop("resize=FALSE is not implemented for class 'sparseSlab'",
         call.=FALSE)
  if(is.null(perm)) return(a)
  stopifnot(length(perm) == 3)
  m <- a$m
  stackdim <- a$stackdim
  newstack <- which(perm == stackdim)
  if(diff(perm[-newstack]) < 0) {
    # transpose indices of sparse matrices
    m <- lapply(m, t)
  }
  sparseSlab(m, stackdim=newstack)
}

"[.sparseSlab" <- local({

  subSlab <- function(x, i,j,k, drop=TRUE, ...) {
    m <- x$m
    stackdim <- x$stackdim
    I <- if(missing(i)) NULL else i
    J <- if(missing(j)) NULL else j
    K <- if(missing(k)) NULL else k
    switch(stackdim,
           {
             ## stackdim = 1
             if(!is.null(I)) m <- listsubset(m, I)
             m <- shootblanks(m, J, K)
           },
           {
             ## stackdim = 2
             if(!is.null(J)) m <- listsubset(m, J)
             m <- shootblanks(m, I, K)
           },
           {
             ## stackdim = 3
             if(!is.null(K)) m <- listsubset(m, K)
             m <- shootblanks(m, I, J)
           })
    if(any(sapply(lapply(m, dim), is.null))) {
      ## This can only happen if a matrix index cbind(i,j) has been used
      y <- matrix(unlist(m), ncol=dim(x)[stackdim])
      dimnames(y)[[2]] <- names(m)
      if(stackdim == 1) y <- t(y)
      return(y[,,drop=drop])
    }
    y <- sparseSlab(m, stackdim=stackdim)

    if(drop) {
      d <- dim(y)
      losing <- (d == 1)
      if(any(losing)) {
        if(losing[stackdim] && !any(losing[-stackdim])) {
          ## result is a sparse matrix
          y <- m[[1]]
        } else {
          ## result is a full matrix or vector
          v <- unlist(lapply(m, as.vector))
          nretained <- sum(!losing)
          if(nretained == 2) {
            ## matrix
            dv <- d[!losing]
            v <- matrix(v, dv[1], dv[2])
            dimnames(v) <- dimnames(y)[!losing]
          } else if(nretained == 1) {
            ## vector of length > 1
            names(v) <- dimnames(y)[!losing][[1]]
          } else {
            ## single value
            v <- unname(v)
          }
          y <- v
        }
      }
    }
    return(y)
  }

  listsubset <- function(x, i) {
    if(!is.vector(i))
      stop("The index for the 'stacking' dimension must not be a matrix",
           call.=FALSE)
    y <- x[i]
    if(any(sapply(y, is.null)))
      stop("subscript out of bounds")
    return(y)
  }

  shootblanks <- function(x, i, j) {
    n <- length(x)
    if(n == 0) return(x)
    y <- vector(mode="list", length=n)
    names(y) <- names(x)
    if(!is.null(i) && !is.null(j)) {
      for(k in seq_len(n))
        y[[k]] <- (x[[k]])[i,j,drop=FALSE]
    } else if(!is.null(i)) {
      ## j is missing or blank
      if(is.vector(i)) {
        for(k in seq_len(n))
          y[[k]] <- (x[[k]])[i,,drop=FALSE]
      } else {
        for(k in seq_len(n))
          y[[k]] <- (x[[k]])[i]
      }
    } else if(!is.null(j)) {
      ## i is missing or blank
      if(is.vector(j)) {
        for(k in seq_len(n))
          y[[k]] <- (x[[k]])[,j,drop=FALSE]
      } else {
        for(k in seq_len(n))
          y[[k]] <- (x[[k]])[j]
      }
    } else return(x)
    return(y)
  }

  subSlab
})

as.array.sparseSlab <- function(x, ...) {
  d <- dim(x)
  dn <- dimnames(x)
  k <- x$stackdim
  # make array with stacks as planes [3rd dim]
  z <- array(unlist(lapply(x$m, as.matrix)),
             dim=c(d[-k], d[k]), dimnames=c(dn[-k], dn[k]))
  # permute
  switch(k,
         { z <- aperm(z, c(3, 1, 2)) },
         { z <- aperm(z, c(1, 3, 2)) },
         { })
  if(!all(dim(z) == d))
    stop("internal error: wrong dimensions")
  return(z)
}

