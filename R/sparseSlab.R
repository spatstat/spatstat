#'
#' sparseSlab.R
#'
#' Sparse 3D arrays or 'slabs' (where one of the dimensions is small)
#' represented as lists of sparse 2D matrices
#' 
#' $Revision: 1.11 $  $Date: 2016/02/22 12:03:55 $
#'

sparseSlab <- function(matlist, stackdim=3) {
  stopifnot(is.list(matlist))
  if(is.null(stackdim)) stackdim <- 3 else stopifnot(stackdim %in% 1:3)
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

emptySparseSlab <- function(typemaker, dim, stackdim=3) {
  stopifnot(length(dim) == 3)
  if(is.null(stackdim)) stackdim <- 3 else stopifnot(stackdim %in% 1:3)
  dmat <- dim[-stackdim]
  matlist <- replicate(dim[stackdim],
                       sparseMatrix(i=integer(0), j=integer(0),
                                    x=typemaker(0), dims=dmat))
  result <- sparseSlab(matlist, stackdim)
  return(result)
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

"[<-.sparseSlab" <- function(x, i, j, k, ..., value) {
  m <- x$m
  stackdim <- x$stackdim
  I <- if(missing(i)) NULL else i
  J <- if(missing(j)) NULL else j
  K <- if(missing(k)) NULL else k
  indices <- list(I,J,K)
  absent <- sapply(indices, is.null)
  involves.stackdim <- !absent[stackdim]
  involves.matrices <- !all(absent[-stackdim])
  dnx <- dimnames(x)
  dimX <- dim(x)
  dimV <- dim(value)
  if(involves.matrices && !involves.stackdim) {
    ## replace subset of each sparse matrix
    for(ii in seq_len(dimX[stackdim])) {
      ind <- indices[-stackdim]
      x$m[[ii]][ind[[1]], ind[[2]]] <- value
    }
  } else if(involves.stackdim && !involves.matrices) {
    ## index selects one or more of the sparse matrices
    stackindex <- positiveIndex(indices[[stackdim]], dnx[[stackdim]])
    nslice <- length(stackindex)
    if(nslice == 1) {
      ## assigning values into a single slice: do it
      x$m[[stackindex]][] <- value
    } else if(nslice > 1) {
      ## assigning values to multiple slices
      if(length(dimV == 3)) {
        if(!all(dimV[-stackdim] == dimX[-stackdim]))
          stop("Replacement dimensions of matrices do not match", call.=FALSE)
        if(dimV[stackdim] == nslice) {
          ## inserting 3D slab into 3D slab
          for(ii in seq_len(nslice)) {
            sii <- stackindex[ii]
            vii <- switch(stackdim,
                          value[ii,,],
                          value[,ii,],
                          value[,,ii])
            if(inherits(vii, "sparseMatrix")) {
              x$m[[ sii ]] <- vii
            } else {
              x$m[[ sii ]] [ ] <- vii[]
            }
          } 
        } else if(dimV[stackdim] == 1) {
          ## replacing slices by the same slice
          value <- switch(stackdim,
                          value[1,,],
                          value[,1,],
                          value[,,1])
          for(ii in seq_len(nslice)) {
            sii <- stackindex[ii]
            if(inherits(value, "sparseMatrix")) {
              x$m[[ sii ]] <- value
            } else {
              x$m[[ sii ]] [ ] <- value[]
            }
          }
        } else 
          stop("I don't know how to match these array indices", call.=FALSE)
      } else if(length(dimV) == 2) {
         ## replacing slices by the same matrix
          for(ii in seq_len(nslice)) 
            x$m[[ stackindex[ii] ]] <- value
      } else if(is.null(dimV) && length(value) == 1) {
        ## replacing slices by a constant 
          for(ii in seq_len(nslice)) 
            x$m[[ stackindex[ii] ]] <- value
      } else stop("Indexing format not understood", call.=FALSE)
    } else stop("Indexing format not understood", call.=FALSE)
  } else
    stop("Sorry, [<-.sparseSlab is not implemented in this case",
         call.=FALSE)
  return(x)
}

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

mapSparseEntries <- function(x, margin, values) {
  # replace the NONZERO entries of sparse matrix or array
  # by values[l] where l is one of the slice indices
  if(inherits(x, "sparseMatrix")) {
    x <- as(x, Class="TsparseMatrix")
    stopifnot(margin %in% 1:2)
    check.nvector(values, dim(x)[margin], things=c("rows","columns")[margin],
                  oneok=TRUE)
    if(length(values) == 1) values <- rep(values, dim(x)[margin])
    i <- x@i + 1L
    j <- x@j + 1L
    yindex <- switch(margin, i, j)
    y <- sparseMatrix(i=i, j=j, x=values[yindex],
                      dims=dim(x), dimnames=dimnames(x))
    return(y)
  }
  stopifnot(inherits(x, "sparseSlab"))
  dims <- dim(x)
  dnx <- dimnames(x)
  stackdim <- x$stackdim
  if(!is.matrix(values)) {
    # vector of values.
    check.nvector(values, dims[margin], oneok=TRUE,
                  things=c("rows", "columns", "planes")[margin])
    if(length(values) == 1) values <- rep(values, dims[margin])
  } else {
    # matrix of values.
    # columns of matrix must match stacking dimension
    # rows of matrix must match 'margin'
    if(margin == stackdim)
      stop("When values are a matrix, margin must not be stacking dimension",
           call.=FALSE)
    if(nrow(values) != dims[margin])
      stop(paste("Number of rows of values", paren(nrow(values)),
                 "does not match array size in margin", paren(dims[margin])),
           call.=FALSE)
    if(ncol(values) != dims[stackdim])
      stop(paste("Number of columns of values", paren(ncol(values)),
                 "does not match array size in stacking dimension",
                 paren(dims[stackdim])),
           call.=FALSE)
  }
  m <- lapply(x$m, as, Class="TsparseMatrix")
  dim.each <- dims[-stackdim]
  dn.each <- dnx[-stackdim]
  for(k in seq_len(dims[stackdim])) {
    mk <- m[[k]]
    i <- mk@i + 1L
    j <- mk@j + 1L
    yindex <- if(margin == stackdim) k else if(margin == 1) i else j
    yvalues <- if(is.matrix(values)) values[yindex, k] else values[yindex]
    m[[k]] <- sparseMatrix(i=i, j=j, x=yvalues,
                           dims=dim.each, dimnames=dn.each)
  }
  x$m <- m
  return(x)
}

anyNA.sparseSlab <- function(x, recursive=FALSE) {
  any(sapply(x$m, anyNA))
}

Ops.sparseSlab <- function(e1,e2=NULL){
  unary <- nargs() == 1L
  if(unary){
    if(!is.element(.Generic, c("!", "-", "+")))
      stop("Unary usage is undefined for this operation for sparse slabs.")
    result <- e1
    result$m <- lapply(result$m, .Generic)
  } else {
    e1val <- if(inherits(e1, "sparseSlab")) e1$m else e1
    e2val <- if(inherits(e2, "sparseSlab")) e2$m else e2
    result <- if(inherits(e1, "sparseSlab")) e1 else e2
    result$m <- mapply(.Generic, e1val, e2val)
  }
  return(result)
}

Math.sparseSlab <- function(x, ...){
  nslice <- length(x$m)
  allstuff <- list(x, ...)
  isslab <- sapply(allstuff, inherits, what="sparseSlab")
  isspam <- sapply(allstuff, inherits, what="sparseMatrix")
  for(i in seq_along(allstuff)) {
    if(isslab[i]) allstuff[[i]] <- allstuff[[i]]$m else
    if(isspam[i]) allstuff[[i]] <- rep(list(allstuff[[i]]), nslice)
  }
  x$m <- do.call(mapply, append(list(.Generic), allstuff))
  return(x)
}

Summary.sparseSlab <- function(..., na.rm) {
  args <- list(...)
  isslab <- sapply(args, inherits, what="sparseSlab")
  mats <-
    if(!any(isslab)) list() else
    Reduce(append, unname(lapply(args[isslab], getElement, name="m")))
  matvals <- lapply(unname(mats), .Generic)
  rslt <- do.call(.Generic, c(matvals, args[!isslab], na.rm=na.rm))
  return(rslt)
}

