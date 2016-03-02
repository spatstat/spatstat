#'
#'    sparsecommon.R
#'
#'  Utilities applicable to both 'sparseSlab' and 'sparse3Darray'
#'
#'  $Revision: 1.1 $  $Date: 2016/03/01 09:02:58 $
#'

#'  .............. completely generic ....................


inside3Darray <- function(d, i, j, k) {
  stopifnot(length(d) == 3)
  if(length(dim(i)) == 2 && missing(j) && missing(k)) {
    stopifnot(ncol(i) == 3)
    j <- i[,2]
    k <- i[,3]
    i <- i[,1]
  }
  ans <- inside.range(i, c(1, d[1])) &
         inside.range(j, c(1, d[2])) &
         inside.range(k, c(1, d[3]))
  return(ans)
}

#'  .............. sparseSlab and sparse3Darray ....................


mapSparseEntries <- function(x, margin, values, conform=TRUE, stackmargin) {
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
  if(inherits(x, "sparse3Darray")) {
    ijk <- cbind(i=x$i, j=x$j, k=x$k)
    if(conform) {
      #' ensure common pattern of sparse values
      #' in dimensions other than 'stackmargin'
      nstack <- dim(x)[stackmargin]
      dup <- duplicated(ijk[,-stackmargin,drop=FALSE])
      ijk <- ijk[!dup, , drop=FALSE]
      npattern <- nrow(ijk)
      #' repeat this pattern in each 'stackmargin' slice
      ijk <- apply(ijk, 2, rep, times=nstack)
      ijk[, stackmargin] <- rep(1:nstack, each=npattern)
    }
    if(!is.matrix(values)) {
      # vector of values matching margin extent
      check.nvector(values, dim(x)[margin],
                    things=c("rows","columns","planes")[margin],
                    oneok=TRUE)
      if(length(values) == 1) values <- rep(values, dim(x)[margin])
      yindex <- ijk[,margin]
      y <- sparse3Darray(i=ijk[,1],
                         j=ijk[,2],
                         k=ijk[,3],
                         x=values[yindex],
                         dims=dim(x), dimnames=dimnames(x))
      return(y)
    } else {
      #' matrix of values.
      force(stackmargin)
      stopifnot(stackmargin != margin) 
      #' rows of matrix must match 'margin'
      if(nrow(values) != dims[margin])
        stop(paste("Number of rows of values", paren(nrow(values)),
                   "does not match array size in margin", paren(dims[margin])),
             call.=FALSE)
      #' columns of matrix must match 'stackmargin'
      if(ncol(values) != dims[stackmargin])
        stop(paste("Number of columns of values", paren(ncol(values)),
                   "does not match array size in stackmargin",
                   paren(dims[stackmargin])),
             call.=FALSE)
      # map
      yindex <- ijk[,margin]
      stackindex <- ijk[,stackmargin]
      y <- sparse3Darray(i=ijk[,1], j=ijk[,2], k=ijk[,3],
                         x=values[cbind(yindex,stackindex)],
                         dims=dim(x), dimnames=dimnames(x))
      return(y)
    } 
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
  m <- x$m
  if(!conform) {
    # use a different pattern of (i,j) for each matrix
    m <- lapply(m, as, Class="TsparseMatrix")
  } else {
    # use a common template pattern of (i,j)
    m <- lapply(m, as, Class="lgCMatrix")
    template <- Reduce("|", m)
    template <- as(template, Class="TsparseMatrix")
    i <- template@i + 1L
    j <- template@j + 1L
  }
  dim.each <- dims[-stackdim]
  dn.each <- dnx[-stackdim]
  for(k in seq_len(dims[stackdim])) {
    if(!conform) {
      mk <- m[[k]]
      i <- mk@i + 1L
      j <- mk@j + 1L
    }
    yindex <- if(margin == stackdim) k else if(margin == 1) i else j
    yvalues <- if(is.matrix(values)) values[yindex, k] else values[yindex]
    m[[k]] <- sparseMatrix(i=i, j=j, x=yvalues,
                           dims=dim.each, dimnames=dn.each)
  }
  x$m <- m
  return(x)
}


applySparseEntries <- local({

  applySparseEntries <- function(x, f, ...) {
    ## apply vectorised function 'f' only to the nonzero entries of 'x'
    if(inherits(x, "sparseMatrix")) {
      x <- applytoxslot(x, f, ...)
    } else if(inherits(x, "sparse3Dmatrix")) {
      x <- applytoxentry(x, f, ...)
    } else if(inherits(x, "sparseSlab")) {
      x$m <- lapply(x$m, applytoxslot, f=f, ...)
    } else {
      x <- f(x, ...)
    }
    return(x)
  }

  applytoxslot <- function(x, f, ...) {
    xx <- x@x
    n <- length(xx)
    xx <- f(xx, ...)
    if(length(xx) != n)
      stop(paste("Function f returned the wrong number of values:",
                 length(xx), "instead of", n),
           call.=FALSE)
    x@x <- xx
    return(x)
  }
  
  applytoxentry <- function(x, f, ...) {
    xx <- x$x
    n <- length(xx)
    xx <- f(xx, ...)
    if(length(xx) != n)
      stop(paste("Function f returned the wrong number of values:",
                 length(xx), "instead of", n),
           call.=FALSE)
    x$x <- xx
    return(x)
  }
  
  applySparseEntries
})

