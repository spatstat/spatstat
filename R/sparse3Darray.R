#'
#' sparse3Darray.R
#'
#' Sparse 3D arrays represented as list(i,j,k,x)
#' 
#' $Revision: 1.33 $  $Date: 2018/05/04 03:14:11 $
#'

sparse3Darray <- function(i=integer(0), j=integer(0), k=integer(0),
                          x=numeric(0),
                          dims=c(max(i),max(j),max(k)),
                          dimnames=NULL, strict=FALSE, nonzero=FALSE) {
  dat <- data.frame(i, j, k, x)
  if(typeof(x) == "complex")
    warn.once("sparse.complex",
              "complex-valued sparse 3D arrays are supported in spatstat,",
              "but complex-valued sparse matrices",
              "are not yet supported by the Matrix package")
  stopifnot(length(dims) == 3)
  dims <- as.integer(dims)
  if(!all(i >= 1 & i <= dims[1])) stop("indices i are outside range")
  if(!all(j >= 1 & j <= dims[2])) stop("indices j are outside range")
  if(!all(k >= 1 & k <= dims[3])) stop("indices k are outside range")
  if(!is.null(dimnames)) {
    stopifnot(is.list(dimnames))
    stopifnot(length(dimnames) == 3)
    notnull <- !sapply(dimnames, is.null)
    dimnames[notnull] <- lapply(dimnames[notnull], as.character)
  }
  if(nonzero || strict) {
    #' drop zeroes
    ok <- (x != RelevantZero(x))
    dat <- dat[ok, , drop=FALSE]
  }
  if(strict) {
    #' arrange in 'R order'
    dat <- dat[with(dat, order(k,j,i)), , drop=FALSE]
    #' duplicates will be adjacent
    dup <- with(dat, c(FALSE, diff(i) == 0 & diff(j) == 0 & diff(k) == 0))
    if(any(dup)) {
      #' accumulate values at the same array location
      retain <- !dup
      newrow <- cumsum(retain)
      newx <- as(tapply(dat$x, newrow, sum), typeof(dat$x))
      newdat <- dat[retain,,drop=FALSE]
      newdat$x <- newx
      dat <- newdat
    }
  }
  result <- append(as.list(dat),
                   list(dim=dims, dimnames=dimnames))
  class(result) <- "sparse3Darray"
  return(result)
}

as.sparse3Darray <- function(x, ...) {
  if(inherits(x, "sparse3Darray")) {
    y <- x
  } else if(inherits(x, c("matrix", "sparseMatrix"))) {
    z <- as(x, Class="TsparseMatrix")
    dn <- dimnames(x)
    dn <- if(is.null(dn)) NULL else c(dn, list(NULL))
    one <- if(length(z@i) > 0) 1L else integer(0)
    y <- sparse3Darray(i=z@i + 1L, j=z@j + 1L, k=one, x=z@x,
                       dims=c(dim(x), 1L), dimnames=dn)
  } else if(is.array(x)) {
    stopifnot(length(dim(x)) == 3)
    dimx <- dim(x)
    if(prod(dimx) == 0) {
      y <- sparse3Darray(, dims=dimx, dimnames=dimnames(x))
    } else {
      ijk <- which(x != RelevantZero(x), arr.ind=TRUE)
      ijk <- cbind(as.data.frame(ijk), x[ijk])
      y <- sparse3Darray(i=ijk[,1L], j=ijk[,2L], k=ijk[,3L], x=ijk[,4L],
                         dims=dimx, dimnames=dimnames(x))
    }
  } else if(inherits(x, "sparseVector")) {
    one <- if(length(x@i) > 0) 1L else integer(0)
    y <- sparse3Darray(i=x@i, j=one, k=one, x=x@x,
                       dims=c(x@length, 1L, 1L))
  } else if(is.null(dim(x)) && is.atomic(x)) {
    n <- length(x)
    dn <- names(x)
    if(!is.null(dn)) dn <- list(dn, NULL, NULL)
    one <- if(n > 0) 1L else integer(0)
    y <- sparse3Darray(i=seq_len(n), j=one, k=one, x=x,
                       dims=c(n, 1L, 1L), dimnames=dn)
    
  } else if(is.list(x) && length(x) > 0) {
    n <- length(x)
    if(all(sapply(x, is.matrix))) {
      z <- Reduce(abind, x)
      y <- as.sparse3Darray(z)
    } else if(all(sapply(x, inherits, what="sparseMatrix"))) {
      dimlist <- unique(lapply(x, dim))
      if(length(dimlist) > 1) stop("Dimensions of matrices do not match")
      dimx <- c(dimlist[[1L]], n)
      dnlist <- lapply(x, dimnames)
      isnul <- sapply(dnlist, is.null)
      dnlist <- unique(dnlist[!isnul])
      if(length(dnlist) > 1) stop("Dimnames of matrices do not match")
      dn <- if(length(dnlist) == 0) NULL else c(dnlist[[1L]], list(NULL))
      for(k in seq_len(n)) {
        mk <- as(x[[k]], "TsparseMatrix")
        kvalue <- if(length(mk@i) > 0) k else integer(0)
        dfk <- data.frame(i=mk@i + 1L, j=mk@j + 1L, k=kvalue, x=mk@x)
        df <- if(k == 1) dfk else rbind(df, dfk)
      }
      y <- sparse3Darray(i=df$i, j=df$j, k=df$k, x=df$x,
                         dims=dimx, dimnames=dn)
    } else {
      warning("I don't know how to convert a list to a sparse array")
      return(NULL)
    }
  } else {          
    warning("I don't know how to convert x to a sparse array")
    return(NULL)
  }
  return(y)
}

dim.sparse3Darray <- function(x) { x$dim }

"dim<-.sparse3Darray" <- function(x, value) {
  stopifnot(length(value) == 3)
  if(!all(inside.range(x$i, c(1, value[1]))))
    stop("indices i are outside new range")
  if(!all(inside.range(x$j, c(1, value[2]))))
    stop("indices j are outside new range")
  if(!all(inside.range(x$k, c(1, value[3]))))
    stop("indices k are outside new range")
  dimx <- dim(x)
  x$dim <- value
  if(!is.null(dimnames(x))) {
    dn <- dimnames(x)
    for(n in 1:3) {
      if(value[n] < dimx[n]) dn[[n]] <- dn[[n]][1:value[n]] else
      if(value[n] > dimx[n]) dn[n] <- list(NULL)
    }
    dimnames(x) <- dn
  }
  return(x)
}

dimnames.sparse3Darray <- function(x) { x$dimnames }

"dimnames<-.sparse3Darray" <- function(x, value) {
  if(!is.list(value)) value <- list(value)
  if(length(value) == 1) value <- rep(value, 3)
  x$dimnames <- value
  return(x)
}

print.sparse3Darray <- function(x, ...) {
  dimx <- dim(x)
  cat("Sparse 3D array of dimensions", paste(dimx, collapse="x"), fill=TRUE)
  if(prod(dimx) == 0)
    return(invisible(NULL))
  dn <- dimnames(x) %orifnull% rep(list(NULL), 3)
  d3 <- dimx[3]
  dn3 <- dn[[3]] %orifnull% as.character(seq_len(d3))
  df <- data.frame(i=x$i, j=x$j, k=x$k, x=x$x)
  pieces <- split(df, factor(df$k, levels=1:d3))
  dim2 <- dimx[1:2]
  dn2 <- dn[1:2]
  if(typeof(x$x) == "complex") {
    splat("\t[Complex-valued sparse matrices are not printable]")
  } else {
    for(k in seq_along(pieces)) {
      cat(paste0("\n\t[ , , ", dn3[k], "]\n\n"))
      Mi <- with(pieces[[k]],
                 sparseMatrix(i=i, j=j, x=x, dims=dim2, dimnames=dn2))
      stuff <- capture.output(eval(Mi))
      #' Remove 'sparse Matrix' header blurb
      stuff <- stuff[-1]
      if(is.blank(stuff[1]))
        stuff <- stuff[-1]
      cat(stuff, sep="\n")
    }
  }
  return(invisible(NULL))
}

aperm.sparse3Darray <- function(a, perm=NULL, resize=TRUE, ...) {
  if(is.null(perm)) return(a)
  stopifnot(length(perm) == 3)
  a <- unclass(a)
  a[c("i", "j", "k")] <- a[c("i", "j", "k")][perm]
  if(resize) {
    a$dim <- a$dim[perm]
    if(length(a$dimnames)==3) a$dimnames <- a$dimnames[perm]
  }
  class(a) <- c("sparse3Darray", class(a))
  return(a)
}

as.array.sparse3Darray <- function(x, ...) {
  zerovalue <- vector(mode=typeof(x$x), length=1L)
  z <- array(zerovalue, dim=dim(x), dimnames=dimnames(x))
  z[cbind(x$i,x$j,x$k)] <- x$x
  return(z)
}

"[.sparse3Darray" <- local({

  Extract <- function(x, i,j,k, drop=TRUE, ...) {
    dimx <- dim(x)
    dn <- dimnames(x) %orifnull% rep(list(NULL), 3)
    if(!missing(i) && length(dim(i)) == 2) {
      ## matrix index
      i <- as.matrix(i)
      if(!(missing(j) && missing(k)))
        stop("If i is a matrix, j and k should not be given", call.=FALSE)
      if(ncol(i) != 3)
        stop("If i is a matrix, it should have 3 columns", call.=FALSE)
      ## start with vector of 'zero' answers of the correct type
      answer <- sparseVector(x=RelevantZero(x$x)[integer(0)],
                             i=integer(0),
                             length=nrow(i))
      ## values outside array return NA
      if(any(bad <- !inside3Darray(dim(x), i)))
        answer[bad] <- NA
      ## if entire array is zero, there is nothing to match
      if(length(x$x) == 0)
        return(answer)
      ## match desired indices to sparse entries
      varies <- (dimx > 1)
      nvary <- sum(varies)
      varying <- which(varies)
      if(nvary == 3) {
        ## ---- older code -----
        ## convert triples of integers to character codes
        #### icode <- apply(i, 1, paste, collapse=",") << is too slow >>
        ## icode <- paste(i[,1], i[,2], i[,3], sep=",")
        ## dcode <- paste(x$i, x$j, x$k, sep=",")
	## ------------------
	m <- matchIntegerDataFrames(i, cbind(x$i, x$j, x$k))
      } else if(nvary == 2) {
        ## effectively a sparse matrix
        ## ---- older code -----
        ## icode <- paste(i[,varying[1]], i[,varying[2]], sep=",")
        ## ijk <- cbind(x$i, x$j, x$k)
        ## dcode <- paste(ijk[,varying[1]], ijk[,varying[2]], sep=",")
	## ------------------
	ijk <- cbind(x$i, x$j, x$k)
	m <- matchIntegerDataFrames(i[,varying,drop=FALSE],
	                            ijk[,varying,drop=FALSE])
      } else if(nvary == 1) {
        ## effectively a sparse vector
        ## ---- older code -----
        ## icode <- i[,varying]
        ## dcode <- switch(varying, x$i, x$j, x$k)
	## ------------------
	m <- match(i[,varying], switch(varying, x$i, x$j, x$k))
      } else {
        ## effectively a single value
        ## ---- older code -----
        ## icode <- rep(1, nrow(i))
        ## dcode <- 1  # since we know length(x$x) > 0
	m <- 1
      }
      ## insert any found elements
      found <- !is.na(m)
      answer[found] <- x$x[m[found]]
      return(answer)
    }
    if(!(missing(i) && missing(j) && missing(k))) {
      I <- grokIndexVector(if(missing(i)) NULL else i, dimx[1], dn[[1]])
      J <- grokIndexVector(if(missing(j)) NULL else j, dimx[2], dn[[2]])
      K <- grokIndexVector(if(missing(k)) NULL else k, dimx[3], dn[[3]])
      IJK <- list(I,J,K)
      if(!all(sapply(lapply(IJK, getElement, name="full"), is.null))) {
        ## indices exceed array bounds; result is a full array containing NA's
        result <- as.array(x)[I$full$i, J$full$j, K$full$k, drop=drop]
        return(result)
      }
      IJK <- lapply(IJK, getElement, name="strict")
      I <- IJK[[1]]
      J <- IJK[[2]]
      K <- IJK[[3]]
      #' number of values to be returned along each margin
      newdims <- sapply(IJK, getElement, name="n")
      #' dimnames of return array
      newdn <- lapply(IJK, getElement, name="s")
      #' find all required data (not necessarily in required order)
      inI <- I$lo
      inJ <- J$lo
      inK <- K$lo
      df <- data.frame(i=x$i, j=x$j, k=x$k, x=x$x)
      use <- with(df, inI[i] & inJ[j] & inK[k])
      df <- df[use, ,drop=FALSE]
      #' contract sub-array to (1:n) * (1:m) * (1:l)
      df <- transform(df,
                      i = cumsum(inI)[i],
                      j = cumsum(inJ)[j],
                      k = cumsum(inK)[k])
      Imap <- I$map
      Jmap <- J$map
      Kmap <- K$map
      if(nrow(df) == 0 || (is.null(Imap) && is.null(Jmap) && is.null(Kmap))) {
        ## return values are already in correct position
        outdf <- df
      } else {
        #' invert map to determine output positions (reorder/repeat entries)
        snI <- seq_len(I$n)
        snJ <- seq_len(J$n)
        snK <- seq_len(K$n)
        imap <- Imap %orifnull% snI
        jmap <- Jmap %orifnull% snJ
        kmap <- Kmap %orifnull% snK
        whichi <- split(seq_along(imap), factor(imap, levels=snI))
        whichj <- split(seq_along(jmap), factor(jmap, levels=snJ))
        whichk <- split(seq_along(kmap), factor(kmap, levels=snK))
        dat.i <- whichi[df$i]
        dat.j <- whichj[df$j]
        dat.k <- whichk[df$k]
        stuff <- mapply(expandwithdata,
                        i=dat.i, j=dat.j, k=dat.k, x=df$x,
                        SIMPLIFY=FALSE)
        outdf <- rbindCompatibleDataFrames(stuff)
      }
      x <- sparse3Darray(i=outdf$i, j=outdf$j, k=outdf$k, x=outdf$x,
                         dims=newdims, dimnames=newdn)
      dimx <- newdims
      dn <- newdn
    }
    if(drop) {
      retain <- (dimx > 1)
      nretain <- sum(retain)
      if(nretain == 2) {
        #' result is a matrix
        retained <- which(retain)
        newi <- getElement(x, name=c("i","j","k")[ retained[1] ])
        newj <- getElement(x, name=c("i","j","k")[ retained[2] ])
        newdim <- dimx[retain]
        newdn <- dn[retain]
        return(sparseMatrix(i=newi, j=newj, x=x$x, dims=newdim, dimnames=newdn))
      } else if(nretain == 1) {
        #' sparse vector
        retained <- which(retain)
        newi <- getElement(x, name=c("i","j","k")[retained])
        #' ensure 'strict' 
        ord <- order(newi)
        newi <- newi[ord]
        newx <- x$x[ord]
        if(any(dup <- c(FALSE, diff(newi) == 0))) {
          retain <- !dup
          ii <- cumsum(retain)
          newi <- newi[retain]
          newx <- as(tapply(newx, ii, sum), typeof(newx))
        }
        x <- sparseVector(x=newx, i=newi, length=dimx[retained])
      } else if(nretain == 0) {
        #' single value
        x <- as.vector(as.array(x))
      }
    }
    return(x)
  }

  expandwithdata <- function(i, j, k, x) {
    z <- expand.grid(i=i, j=j, k=k)
    if(nrow(z) > 0)
      z$x <- x
    return(z)
  }

  Extract
})


rbindCompatibleDataFrames <- function(x) {
  #' faster version of Reduce(rbind, x) when entries are known to be compatible
  nama2 <- colnames(x[[1]])
  y <- vector(mode="list", length=length(nama2))
  names(y) <- nama2
  for(nam in nama2)
    y[[nam]] <- unlist(lapply(x, getElement, name=nam))
  return(as.data.frame(y))
}


"[<-.sparse3Darray" <- function(x, i, j, k, ..., value) {
  dimx <- dim(x)
  dn <- dimnames(x) %orifnull% rep(list(NULL), 3)
  #' interpret indices
  if(!missing(i) && length(dim(i)) == 2) {
    ## matrix index
    ijk <- as.matrix(i)
    if(!(missing(j) && missing(k)))
      stop("If i is a matrix, j and k should not be given", call.=FALSE)
    if(ncol(ijk) != 3)
      stop("If i is a matrix, it should have 3 columns", call.=FALSE)
    if(!all(inside3Darray(dimx, i)))
      stop("Some indices lie outside array limits", call.=FALSE)
    if(nrow(ijk) == 0)
      return(x) # no items to replace
    ## assemble data frame
    xdata <- data.frame(i=x$i, j=x$j, k=x$k, x=x$x)
    ## match xdata into ijk (not necessarily the first match in original order)
    m <- matchIntegerDataFrames(xdata[,1:3,drop=FALSE], ijk)
    ## ------- OLDER VERSION: --------
    ## convert triples of integers to character codes
    ## icode <- apply(ijk, 1, paste, collapse=",") << is too slow >>
    ## icode <- paste(ijk[,1], ijk[,2], ijk[,3], sep=",")
    ## xcode <- paste(x$i, x$j, x$k, sep=",")
    ##  m <- match(xcode, icode)
    ## -------------------------------
    ## remove any matches, retaining only data that do not match 'i'
    xdata <- xdata[is.na(m), , drop=FALSE]   # sic
    ## ensure replacement value is vector-like
    value <- as.vector(value)
    nv <- length(value)
    if(nv != nrow(i) && nv != 1)
      stop(paste("Number of items to replace", paren(nrow(i)),
                 "does not match number of items given", paren(nv)),
           call.=FALSE)
    vdata <- data.frame(i=ijk[,1], j=ijk[,2], k=ijk[,3], x=value)
    ## combine
    ydata <- rbind(xdata, vdata)
    y <- with(ydata, sparse3Darray(i=i,j=j,k=k,x=x,
                                   dims=dimx, dimnames=dn, strict=TRUE))
    return(y)
  }
  I <- grokIndexVector(if(missing(i)) NULL else i, dimx[1], dn[[1]])
  J <- grokIndexVector(if(missing(j)) NULL else j, dimx[2], dn[[2]])
  K <- grokIndexVector(if(missing(k)) NULL else k, dimx[3], dn[[3]])
  IJK <- list(I,J,K)
  if(!all(sapply(lapply(IJK, getElement, name="full"), is.null))) {
    warning("indices exceed array bounds; using full array", call.=FALSE)
    x <- as.array(x)
    x[I$full$i, J$full$j, K$full$k] <- value
    x <- as.sparse3Darray(x)
    return(x)
  }
  IJK <- lapply(IJK, getElement, name="strict")
  if(all(sapply(IJK, getElement, name="nind") == 0)) {
    # no elements are indexed
    return(x)
  }
  I <- IJK[[1]]
  J <- IJK[[2]]
  K <- IJK[[3]]
  #' extract current array entries
  xdata <- data.frame(i=x$i, j=x$j, k=x$k, x=x$x)
  #' identify data volume that will be overwritten
  inI <- I$lo
  inJ <- J$lo
  inK <- K$lo
  #' remove data that will be overwritten
  retain <- !with(xdata, inI[i] & inJ[j] & inK[k])
  xdata <- xdata[retain,,drop=FALSE]
  #' expected dimensions of 'value' implied by indices
  dimVshould <- sapply(IJK, getElement, name="nind")
  dimV <- dim(value)
  if(length(dimV) == 3) {
    #' both source and destination are 3D
    if(all(dimVshould == dimV)) {
      #' replace 3D block by 3D block of same dimensions
      value <- as.sparse3Darray(value)
      vdata <- data.frame(i=value$i, j=value$j, k=value$k, x=value$x)
      # determine positions of replacement data in original array
      vdata <- transform(vdata,
                         i=replacementIndex(i, I),
                         j=replacementIndex(j, J),
                         k=replacementIndex(k, K))
    } else
      stop(paste("Replacement value has wrong dimensions:",
                 paste(dimV, collapse="x"),
                 "instead of",
                 paste(dimVshould, collapse="x")),
           call.=FALSE)
  } else if(is.null(dimV)) {
    #' replacement value is a vector or sparseVector
    value <- as(value, "sparseVector")
    iv <- value@i
    xv <- value@x
    nv <- value@length
    collapsing <- (dimVshould == 1)
    realdim <- sum(!collapsing)
    if(nv == 1) {
      #' replacement value is a constant
      value <- as.vector(value[1])
      if(identical(value, RelevantZero(x$x))) {
        #' assignment causes relevant entries to be set to zero;
        #' these entries have already been deleted from 'xdata';
        #' nothing to add
        vdata <- data.frame(i=integer(0), j=integer(0), k=integer(0),
                            x=x$x[integer(0)])
      } else {
        #' replicate the constant
        vdata <- expand.grid(i=I$i, j=J$i, k=K$i, x=as.vector(value[1]))
      }
    } else if(realdim == 0) {
        stop(paste("Replacement value has too many entries:",
                   nv, "instead of 1"),
             call.=FALSE)
    } else if(realdim == 1) {
      theindex <- which(!collapsing)
      # target slice is one-dimensional
      if(nv != dimVshould[theindex]) 
        stop(paste("Replacement value has wrong number of entries:",
                   nv, "instead of", dimVshould[theindex]),
             call.=FALSE)
      newpos <- replacementIndex(iv, IJK[[theindex]])
      vdata <- switch(theindex,
                      data.frame(i=newpos, j=J$i,    k=K$i,     x=xv),
                      data.frame(i=I$i,    j=newpos, k=K$i,     x=xv),
                      data.frame(i=I$i,    j=J$i,    k=newpos,  x=xv))
    } else {
      # target slice is two-dimensional
      sdim <- dimVshould[!collapsing]
      sd1 <- sdim[1]
      sd2 <- sdim[2]
      if(nv != sd1)
        stop(paste("Length of replacement vector", paren(nv),
                   "does not match dimensions of array subset",
                   paren(paste(dimVshould, collapse="x"))),
             call.=FALSE)
      firstindex <- which(!collapsing)[1]
      secondindex <- which(!collapsing)[2]
      pos1 <- replacementIndex(iv, IJK[[firstindex]])
      pos2 <- replacementIndex(seq_len(sd2), IJK[[secondindex]])
      xv   <- rep(xv, sd2)
      pos1 <- rep(pos1, sd2)
      pos2 <- rep(pos2, each=length(pos1))
      pos3 <- if(length(pos1)) IJK[[which(collapsing)]]$i else integer(0)
      vdata <- data.frame(i=pos3, j=pos3, k=pos3, x=xv)
      vdata[,firstindex] <- pos1
      vdata[,secondindex] <- pos2
    }
  } else if(identical(dimVshould[dimVshould > 1],  dimV[dimV > 1])) {
    #' lower dimensional sets of the same dimension
    value <- value[drop=TRUE]
    dimV <- dim(value)
    dropping <- (dimVshould == 1)
    if(length(dimV) == 2) {
      value <- as(value, "TsparseMatrix")
      iv <- value@i + 1L
      jv <- value@j + 1L
      xv <- value@x
      firstindex <- which(!dropping)[1]
      secondindex <- which(!dropping)[2]
      pos1 <- replacementIndex(iv, IJK[[firstindex]])
      pos2 <- replacementIndex(jv, IJK[[secondindex]])
      pos3 <- if(length(pos1)) IJK[[which(dropping)]]$i else integer(0)
      vdata <- data.frame(i=pos3, j=pos3, k=pos3, x=xv)
      vdata[,firstindex] <- pos1
      vdata[,secondindex] <- pos2
    } else {
      value <- as(value, "sparseVector")
      iv <- value@i
      xv <- value@x
      vdata <- data.frame(i=if(dropping[1]) I$i else replacementIndex(iv, I),
                          j=if(dropping[2]) J$i else replacementIndex(iv, J),
                          k=if(dropping[3]) K$i else replacementIndex(iv, K),
                          x=xv)
    }
  } else
    stop(paste("Replacement value has wrong dimensions:",
               paste(dimV, collapse="x"),
               "instead of",
               paste(dimVshould, collapse="x")),
         call.=FALSE)
    
  ## combine
  if(nrow(vdata) > 0)
    xdata <- rbind(xdata, vdata)
  y <- with(xdata, sparse3Darray(i=i,j=j,k=k,x=x,
                                 dims=dimx, dimnames=dn, strict=TRUE))
  return(y)
}

bind.sparse3Darray <- function(A,B,along) {
  A <- as.sparse3Darray(A)
  B <- as.sparse3Darray(B)
  check.1.integer(along)
  stopifnot(along %in% 1:3)
  dimA <- dim(A)
  dimB <- dim(B)
  if(!all(dimA[-along] == dimB[-along]))
    stop("dimensions of A and B do not match")
  dimC <- dimA
  dimC[along] <- dimA[along] + dimB[along]
  # extract data
  Adf <- SparseEntries(A)
  Bdf <- SparseEntries(B)
  # realign 'B' coordinate
  Bdf[,along] <- Bdf[,along] + dimA[along]
  # combine
  C <- EntriesToSparse(rbind(Adf, Bdf), dimC)
  # add dimnames
  dnA <- dimnames(A)
  dnB <- dimnames(B)
  if(!is.null(dnA) || !is.null(dnB)) {
    if(length(dnA) != 3) dnA <- rep(list(NULL), 3)
    if(length(dnB) != 3) dnB <- rep(list(NULL), 3)
    dnC <- dnA
    dnC[[along]] <- c(dnA[[along]] %orifnull% rep("", dimA[along]),
                      dnB[[along]] %orifnull% rep("", dimB[along]))
    dimnames(C) <- dnC
  }
  return(C)
}


anyNA.sparse3Darray <- function(x, recursive=FALSE) {
  anyNA(x$x)
}

RelevantZero <- function(x) vector(mode=typeof(x), length=1L)
isRelevantZero <- function(x) identical(x, RelevantZero(x))
RelevantEmpty <- function(x) vector(mode=typeof(x), length=0L)

unionOfSparseIndices <- function(A, B) {
  #' A, B are data frames of indices i, j, k
  ijk <- unique(rbind(A, B))
  colnames(ijk) <- c("i", "j", "k")
  return(ijk)
}
  
Ops.sparse3Darray <- function(e1,e2=NULL){
  if(nargs() == 1L) {
    switch(.Generic,
           "!" = {
             result <- do.call(.Generic, list(as.array(e1)))
           },
           "-" = ,
           "+" = {
             result <- e1
             result$x <- do.call(.Generic, list(e1$x))
           },
           stop(paste("Unary", sQuote(.Generic),
                      "is undefined for sparse 3D arrays."), call.=FALSE))
    return(result)
  }
  # binary operation
  # Decide whether full or sparse
  elist <- list(e1, e2)
  isfull <- sapply(elist, inherits, what=c("matrix", "array"))
  if(any(isfull) &&
     any(sapply(lapply(elist[isfull], dim), prod) > 1)) {
    # full array
    n1 <- length(dim(e1))
    n2 <- length(dim(e2))
    e1 <- if(n1 == 3) as.array(e1) else
          if(n1 == 2) as.matrix(e1) else as.vector(as.matrix(as.array(e1)))
    e2 <- if(n2 == 3) as.array(e2) else
          if(n2 == 2) as.matrix(e2) else as.vector(as.matrix(as.array(e2)))
    result <- do.call(.Generic, list(e1, e2))
    return(result)
  }
  # sparse result (usually)
  e1 <- as.sparse3Darray(e1)
  e2 <- as.sparse3Darray(e2)
  dim1 <- dim(e1)
  dim2 <- dim(e2)
  mode1 <- typeof(e1$x)
  mode2 <- typeof(e2$x)
  zero1 <- vector(mode=mode1, length=1L)
  zero2 <- vector(mode=mode2, length=1L)
  
  if(prod(dim1) == 1) {
    ## e1 is constant
    e1 <- as.vector(as.array(e1))
    z12 <- do.call(.Generic, list(e1, zero2))
    if(!isRelevantZero(z12)) {
      # full matrix/array will be generated
      result <- do.call(.Generic, list(e1, as.array(e2)[drop=TRUE]))
    } else {
      # sparse 
      result <- e2
      result$x <- do.call(.Generic, list(e1, e2$x))
    }
    return(result)
  }

  if(prod(dim2) == 1) {
    ## e2 is constant
    e2 <- as.vector(as.array(e2))
    z12 <- do.call(.Generic, list(zero1, e2))
    if(!isRelevantZero(z12)) {
      # full matrix/array will be generated
      result <- do.call(.Generic, list(as.array(e1)[drop=TRUE], e2))
    } else {
      # sparse 
      result <- e1
      result$x <- do.call(.Generic, list(e1$x, e2))
    }
    return(result)
  }
  
  z12 <- do.call(.Generic, list(zero1, zero2))
  if(!isRelevantZero(z12)) {
    #' Result is an array
    e1 <- as.array(e1)
    e2 <- as.array(e2)
    result <- do.call(.Generic, list(e1, e2))
    return(result)
  }

  # Result is sparse
  if(identical(dim1, dim2)) {
    #' extents are identical
    ijk1 <- SparseIndices(e1)
    ijk2 <- SparseIndices(e2)
    if(identical(ijk1, ijk2)) {
      #' patterns of nonzero entries are identical
      ijk <- ijk1
      values <- do.call(.Generic, list(e1$x, e2$x))
    } else {			   
      #' different patterns of nonzero entries
      ijk <- unionOfSparseIndices(ijk1, ijk2)
      values <- as.vector(do.call(.Generic, list(e1[ijk], e2[ijk])))
    }			      
    dn <- dimnames(e1) %orifnull% dimnames(e2)
    result <- sparse3Darray(i=ijk$i, j=ijk$j, k=ijk$k, x=values,
                              dims=dim1, dimnames=dn, strict=TRUE)
    return(result)
  }

  drop1 <- (dim1 == 1)
  drop2 <- (dim2 == 1)
  if(!any(drop1 & !drop2) && identical(dim1[!drop2], dim2[!drop2])) {
    #' dim2 is a slice of dim1
    ijk1 <- data.frame(i=e1$i, j=e1$j, k=e1$k)
    ijk2 <- data.frame(i=e2$i, j=e2$j, k=e2$k)
    expanding <- which(drop2 & !drop1)
    if(length(expanding) == 1) {
      n <- dim1[expanding]
      m <- nrow(ijk2)
      ijk2 <- as.data.frame(lapply(ijk2, rep, times=n))
      ijk2[,expanding] <- rep(seq_len(n), each=m)
      ijk <- unionOfSparseIndices(ijk1, ijk2)
      ijkdrop <- ijk
      if(nrow(ijkdrop) > 0) ijkdrop[,expanding] <- 1
      xout <- do.call(.Generic, list(e1[ijk], e2[ijkdrop]))
      result <- sparse3Darray(i=ijk[,1L], j=ijk[,2L], k=ijk[,3L],
                              x=as.vector(xout),
                              dims=dim1, dimnames=dimnames(e1), strict=TRUE)
      return(result)
    }
  }

  if(!any(drop2 & !drop1) && identical(dim2[!drop1], dim1[!drop1])) {
    #' dim1 is a slice of dim2
    ijk1 <- data.frame(i=e1$i, j=e1$j, k=e1$k)
    ijk2 <- data.frame(i=e2$i, j=e2$j, k=e2$k)
    expanding <- which(drop1 & !drop2)
    if(length(expanding) == 1) {
      n <- dim2[expanding]
      m <- nrow(ijk1)
      ijk1 <- as.data.frame(lapply(ijk1, rep, times=n))
      ijk1[,expanding] <- rep(seq_len(n), each=m)
      ijk <- unionOfSparseIndices(ijk1, ijk2)
      ijkdrop <- ijk
      if(nrow(ijkdrop) > 0) ijkdrop[,expanding] <- 1L
      xout <- do.call(.Generic, list(e1[ijkdrop], e2[ijk]))
      result <- sparse3Darray(i=ijk[,1L], j=ijk[,2L], k=ijk[,3L],
                              x=as.vector(xout),
                              dims=dim2, dimnames=dimnames(e2), strict=TRUE)
      return(result)
    }
  }

  if(all(drop1[-1]) && dim1[1L] == dim2[1L]) {
    #' e1 is a (sparse) vector matching the first extent of e2
    if(.Generic %in% c("*", "&")) {
      # result is sparse
      ijk <- data.frame(i=e2$i, j=e2$j, k=e2$k)
      ones <- rep(1L, nrow(ijk))
      i11 <- data.frame(i=e2$i, j=ones, k=ones)
      xout <- do.call(.Generic, list(e1[i11], e2[ijk]))
      result <- sparse3Darray(i=ijk[,1L], j=ijk[,2L], k=ijk[,3L],
                              x=as.vector(xout),
                              dims=dim2, dimnames=dimnames(e2), strict=TRUE)
    } else {
      # result is full array
      e1 <- as.array(e1)[,,,drop=TRUE]
      e2 <- as.array(e2)
      result <- do.call(.Generic, list(e1, e2))
    }
    return(result)
  }
  
  stop(paste("Non-conformable arrays:",
             paste(dim1, collapse="x"), "and", paste(dim2, collapse="x")),
       call.=FALSE)
}

Math.sparse3Darray <- function(x, ...){
  z <- RelevantZero(x$x)
  fz <- do.call(.Generic, list(z))
  if(!isRelevantZero(fz)) {
    # result is a full array
    result <- do.call(.Generic, list(as.array(x), ...))
    return(result)
  }
  x$x <- do.call(.Generic, list(x$x))
  return(x)
}

Summary.sparse3Darray <- function(..., na.rm=FALSE) {
  argh <- list(...)
  is3D <- sapply(argh, inherits, what="sparse3Darray")
  if(any(is3D)) {
    xvalues   <- lapply(argh[is3D], getElement, name="x")
    fullsizes <- sapply(lapply(argh[is3D], dim), prod)
    argh[is3D] <- xvalues
    #' zero entry should be appended if and only if there are any empty cells
    zeroes <- lapply(xvalues, RelevantZero)
    zeroes <- zeroes[lengths(xvalues) < fullsizes]
    argh <- append(argh, zeroes)
  }
  rslt <- do.call(.Generic, append(argh, list(na.rm=na.rm)))
  return(rslt)
}


SparseIndices <- function(x) {
  #' extract indices of entries of sparse vector/matrix/array
  nd <- length(dim(x))
  if(nd > 3)
    stop("Arrays of more than 3 dimensions are not supported", call.=FALSE)
  if(nd == 0 || nd == 1) {
    x <- as(x, "sparseVector")
    df <- data.frame(i=x@i)
  } else if(nd == 2) {
    x <- as(x, "TsparseMatrix")
    df <- data.frame(i=x@i + 1L, j=x@j + 1L)
  } else if(nd == 3) {
    x <- as.sparse3Darray(x)
    df <- data.frame(i=x$i, j=x$j, k=x$k)
  }
  return(df)
}

SparseEntries <- function(x) {
  #' extract entries of sparse vector/matrix/array
  nd <- length(dim(x))
  if(nd > 3)
    stop("Arrays of more than 3 dimensions are not supported", call.=FALSE)
  if(nd == 0 || nd == 1) {
    x <- as(x, "sparseVector")
    df <- data.frame(i=x@i, x=x@x)
  } else if(nd == 2) {
    x <- as(x, "TsparseMatrix")
    df <- data.frame(i=x@i + 1L, j=x@j + 1L, x=x@x)
  } else if(nd == 3) {
    x <- as.sparse3Darray(x)
    df <- data.frame(i=x$i, j=x$j, k=x$k, x=x$x)
  }
  return(df)
}

EntriesToSparse <- function(df, dims) {
  #' convert data frame of indices and values
  #' to sparse vector/matrix/array
  nd <- length(dims)
  if(nd == 0)
    return(with(df, as(sum(x), typeof(x))))
  sn <- seq_len(nd)
  colnames(df)[sn] <- c("i","j","k")[sn]
  if(nd == 1) {
    #' sparse vector: duplicate entries not allowed
    df <- df[with(df, order(i)), , drop=FALSE]
    dup <- c(FALSE, with(df, diff(i) == 0))
    if(any(dup)) {
      #' accumulate values at the same array location
      first <- !dup
      newi <- cumsum(first)
      newx <- as(tapply(df$x, newi, sum), typeof(df$x))
      df <- data.frame(i=newi[first], x=newx)
    }
    result <- with(df, sparseVector(i=i, x=x, length=dims))
  } else if(nd == 2) {
    result <- with(df, sparseMatrix(i=i, j=j, x=x, dims=dims))
  } else if(nd == 3) {
    result <- with(df, sparse3Darray(i=i, j=j, k=k, x=x, dims=dims))
  }
  return(result)
}

evalSparse3Dentrywise <- function(expr, envir) {
  ## DANGER: this assumes all sparse arrays in the expression
  ##         have the same pattern of nonzero elements!
  e <- as.expression(substitute(expr))
  ## get names of all variables in the expression
  varnames <- all.vars(e)
  allnames <- all.names(e, unique=TRUE)
  funnames <- allnames[!(allnames %in% varnames)]
  if(length(varnames) == 0)
    stop("No variables in this expression")
  ## get the values of the variables
  if(missing(envir)) {
    envir <- parent.frame() # WAS: sys.parent()
  } else if(is.list(envir)) {
    envir <- list2env(envir, parent=parent.frame())
  }
  vars <- mget(varnames, envir=envir, inherits=TRUE, ifnotfound=list(NULL))
  funs <- mget(funnames, envir=envir, inherits=TRUE, ifnotfound=list(NULL))
  ## find out which variables are sparse3Darray
  isSpud <- sapply(vars, inherits, what="sparse3Darray")
  if(!any(isSpud))
    stop("No sparse 3D arrays in this expression")
  spuds <- vars[isSpud]
  template <- spuds[[1L]]
  ## replace each array by its entries, and evaluate
  spudvalues <- lapply(spuds, getElement, name="x")
  ## minimal safety check
  if(length(unique(lengths(spudvalues))) > 1)
    stop("Different numbers of sparse entries", call.=FALSE)
  vars[isSpud] <- spudvalues
  v <- eval(e, append(vars, funs))
  ## reshape as 3D array
  result <- sparse3Darray(x=v,
  	                  i=template$i,
  	                  j=template$j,
  	                  k=template$k,
			  dims=dim(template),
			  dimnames=dimnames(template))
  return(result)
}
