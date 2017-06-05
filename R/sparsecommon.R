#'
#'    sparsecommon.R
#'
#'  Utilities for sparse arrays
#'
#'  $Revision: 1.6 $  $Date: 2017/06/05 10:31:58 $
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

#'  .............. depends on Matrix package ................

sparseVectorCumul <- function(x, i, length) {
  #' extension of 'sparseVector' to allow repeated indices 
  #'   (the corresponding entries are added)
  z <- tapply(x, list(factor(i, levels=1:length)), sum)
  z <- z[!is.na(z)]
  sparseVector(i=as.integer(names(z)), x=as.numeric(z), length=length)
}

#'  .............. code that mentions sparse3Darray ................

expandSparse <- function(x, n, across) {
  #' x is a sparse vector/matrix; replicate it 'n' times
  #' and form a sparse matrix/array
  #' in which each slice along the 'across' dimension is identical to 'x'
  #' Default is across = length(dim(x)) + 1
  check.1.integer(n)
  stopifnot(n >= 1)
  dimx <- dim(x)
  if(is.null(dimx)) {
    if(inherits(x, "sparseVector")) dimx <- x@length else
    if(is.vector(x)) dimx <- length(x) else
    stop("Format of x is not understood", call.=FALSE)
  }
  nd <- length(dimx)
  if(missing(across)) across <- nd + 1L else {
    check.1.integer(across)
    if(!(across %in% (1:(nd+1L))))
      stop(paste("Argument 'across' must be an integer from 1 to", nd+1L),
           call.=FALSE)
  }
  if(nd == 1) {
    if(inherits(x, "sparseVector")) {
      m <- length(x@x)
      y <- switch(across,
                  sparseMatrix(i=rep(1:n, times=m),
		               j=rep(x@i, each=n),
			       x=rep(x@x, each=n),
			       dims=c(n, dimx)),
                  sparseMatrix(i=rep(x@i, each=n),
		  	       j=rep(1:n, times=m),
			       x=rep(x@x, each=n),
			       dims=c(dimx, n)))
    } else {
      y <- switch(across,
                  outer(1:n, x, function(a,b) b),
                  outer(x, 1:n, function(a,b) a))
    }
  } else if(nd == 2) {
    if(inherits(x, "sparseMatrix")) {
      z <- as(x, "TsparseMatrix")
      m <- length(z@x)
      y <- switch(across,
                  sparse3Darray(i=rep(1:n, times=m),
		                j=rep(z@i + 1L, each=n),
				k=rep(z@j + 1L, each=n),
				x=rep(z@x, each=n),
				dims=c(n, dimx)),
                  sparse3Darray(i=rep(z@i + 1L, each=n),
		                j=rep(1:n, times=m),
				k=rep(z@j + 1L, each=n),
				x=rep(z@x, each=n),
				dims=c(dimx[1], n, dimx[2])),
                  sparse3Darray(i=rep(z@i + 1L, each=n),
		                j=rep(z@j + 1L, each=n),
				k=rep(1:n, times=m),
				x=rep(z@x, each=n),
				dims=c(dimx, n)))
    } else stop("Not yet implemented for full arrays")
  } else 
     stop("Not implemented for arrays of more than 2 dimensions", call.=FALSE)
  return(y)
}

mapSparseEntries <- function(x, margin, values, conform=TRUE, across) {
  # replace the NONZERO entries of sparse matrix or array
  # by values[l] where l is one of the slice indices
  dimx <- dim(x)
  if(is.null(dimx)) {
    if(inherits(x, "sparseVector")) dimx <- x@length else
    if(is.vector(x)) dimx <- length(x) else
    stop("Format of x is not understood", call.=FALSE)
  }
  if(inherits(x, "sparseMatrix")) {
    x <- as(x, Class="TsparseMatrix")
    if(length(x$i) == 0) {
      # no entries
      return(x)
    }
    check.1.integer(margin)
    stopifnot(margin %in% 1:2)
    check.nvector(values, dimx[margin], things=c("rows","columns")[margin],
                  oneok=TRUE)
    if(length(values) == 1) values <- rep(values, dimx[margin])
    i <- x@i + 1L
    j <- x@j + 1L
    yindex <- switch(margin, i, j)
    y <- sparseMatrix(i=i, j=j, x=values[yindex],
                      dims=dimx, dimnames=dimnames(x))
    return(y)
  }
  if(inherits(x, "sparse3Darray")) {
    if(length(x$i) == 0) {
      # no entries
      return(x)
    }
    ijk <- cbind(i=x$i, j=x$j, k=x$k)
    if(conform) {
      #' ensure common pattern of sparse values
      #' in each slice on 'across' margin
      nslice <- dimx[across]
      #' pick one representative of each equivalence class
      ## ---- old code ---------
      ## dup <- duplicated(ijk[,-across,drop=FALSE])
      ## ijk <- ijk[!dup, , drop=FALSE]
      ## ---------------------
      use <- representativeRows(ijk[,-across,drop=FALSE])
      ijk <- ijk[use, , drop=FALSE]
      ## 
      npattern <- nrow(ijk)
      #' repeat this pattern in each 'across' slice
      ijk <- apply(ijk, 2, rep, times=nslice)
      ijk[, across] <- rep(seq_len(nslice), each=npattern)
    }
    if(is.vector(values)) {
      # vector of values matching margin extent
      check.nvector(values, dimx[margin],
                    things=c("rows","columns","planes")[margin],
                    oneok=TRUE)
      if(length(values) == 1) values <- rep(values, dimx[margin])
      yindex <- ijk[,margin]
      y <- sparse3Darray(i=ijk[,1],
                         j=ijk[,2],
                         k=ijk[,3],
                         x=values[yindex],
                         dims=dimx, dimnames=dimnames(x))
      return(y)
    } else if(is.matrix(values)) {
      #' matrix of values.
      force(across)
      stopifnot(across != margin) 
      #' rows of matrix must match 'margin'
      if(nrow(values) != dimx[margin])
        stop(paste("Number of rows of values", paren(nrow(values)),
                   "does not match array size in margin", paren(dimx[margin])),
             call.=FALSE)
      #' columns of matrix must match 'across'
      if(ncol(values) != dimx[across])
        stop(paste("Number of columns of values", paren(ncol(values)),
                   "does not match array size in 'across'",
                   paren(dimx[across])),
             call.=FALSE)
      # map
      yindex <- ijk[,margin]
      zindex <- ijk[,across]
      y <- sparse3Darray(i=ijk[,1], j=ijk[,2], k=ijk[,3],
                         x=values[cbind(yindex,zindex)],
                         dims=dimx, dimnames=dimnames(x))
      return(y)
    } else stop("Format of values not understood", call.=FALSE)
  }
  stop("Format of x not understood", call.=FALSE)
}


applySparseEntries <- local({

  applySparseEntries <- function(x, f, ...) {
    ## apply vectorised function 'f' only to the nonzero entries of 'x'
    if(inherits(x, "sparseMatrix")) {
      x <- applytoxslot(x, f, ...)
    } else if(inherits(x, "sparse3Dmatrix")) {
      x <- applytoxentry(x, f, ...)
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

