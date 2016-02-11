##
## hypersub.R
##
##
##  subset operations for hyperframes
##
##  $Revision: 1.24 $    $Date: 2016/02/11 10:17:12 $
##

"[.hyperframe" <- function(x, i, j, drop=FALSE, strip=drop, ...) {
  x <- unclass(x)
  if(!missing(i)) {
    y <- x
    y$df     <- x$df[i, , drop=FALSE]
    y$ncases <- nrow(y$df)
    y$hypercolumns <- lapply(x$hypercolumns, "[", i=i)
    x <- y
  }
  if(!missing(j)) {
    y <- x
    patsy <- seq_len(y$nvars)
    names(patsy) <- y$vname
    jj <- patsy[j]
    names(jj) <- NULL
    y$nvars <- length(jj)
    y$vname <- vname <- x$vname[jj]
    y$vtype <- vtype <- x$vtype[jj]
    y$vclass <- x$vclass[jj]
    if(ncol(x$df) != 0) 
      y$df    <- x$df[ , vname[vtype == "dfcolumn"], drop=FALSE]
    y$hyperatoms <- x$hyperatoms[ vname[ vtype == "hyperatom" ]]
    y$hypercolumns <- x$hypercolumns[ vname [ vtype == "hypercolumn" ] ]
    x <- y
  }
  if(drop) {
    nrows <- x$ncases
    ncols <- x$nvars
    if(nrows == 1 && ncols == 1 && strip) {
      ## return a single object 
      y <- switch(as.character(x$vtype),
                  dfcolumn    = x$df[, , drop=TRUE],
                  hypercolumn = (x$hypercolumns[[1]])[[1]],
                  hyperatom   = x$hyperatoms[[1]])
      return(y)
    } else if(nrows == 1) {
      ## return the row as a vector or a list
      if(strip && all(x$vtype == "dfcolumn"))
        return(x$df[ , , drop=TRUE])
      n <- x$nvars
      y <- vector(mode="list", length=n)
      names(y) <- nama <- x$vname
      for(i in seq_len(n)) {
        nami <- nama[i]
        y[[i]] <- switch(as.character(x$vtype[i]),
                         dfcolumn = x$df[ , nami, drop=TRUE],
                         hyperatom = x$hyperatoms[[nami]],
                         hypercolumn = (x$hypercolumns[[nami]])[[1]]
                         )
      }
      return(as.solist(y, demote=TRUE))
    } else if(ncols == 1) {
      ## return a column as an 'anylist'/'solist' or a vector
      switch(as.character(x$vtype),
             dfcolumn = {
               return(x$df[, , drop=TRUE])
             },
             hypercolumn = {
               y <- as.solist(x$hypercolumns[[1]], demote=TRUE)
               names(y) <- row.names(x$df)
               return(y)
             },
             hyperatom = {
               ## replicate it to make a hypercolumn
               ha <- x$hyperatoms[1]
               names(ha) <- NULL
               hc <- rep.int(ha, x$ncases)
               hc <- as.solist(hc, demote=TRUE)
               names(hc) <- row.names(x$df)
               return(hc)
             }
           )
    }
  }
  class(x) <- c("hyperframe", class(x))
  return(x)
}

"$.hyperframe" <- function(x,name) {
  m <- match(name, unclass(x)$vname)
  if(is.na(m))
    return(NULL)
  return(x[, name, drop=TRUE, strip=FALSE])
}

"$<-.hyperframe" <- function(x, name, value) {
  y <- as.list(x)
  if(is.hyperframe(value)) {
    if(ncol(value) == 1) {
      y[name] <- as.list(value)
    } else {
      y <- insertinlist(y, name, as.list(value))
    }
  } else {
    dfcol <- is.atomic(value) && (is.vector(value) || is.factor(value))
    if(!dfcol && !is.null(value))
      value <- as.list(value)
    y[[name]] <- value
  }
  z <- do.call(hyperframe, append(y, list(row.names=row.names(x),
                                            stringsAsFactors=FALSE)))
  return(z)
}

"[<-.hyperframe" <- 
function (x, i, j, value)
{
  sumry <- summary(x)
  colnam <- sumry$col.names
  dimx <- sumry$dim
  igiven <- !missing(i)
  jgiven <- !missing(j)
  if(!igiven) i <- seq_len(dimx[1])
  if(!jgiven) j <- seq_len(dimx[2])
#  singlerow    <- ((is.integer(i) && length(i) == 1 && i > 0)
#                   || (is.character(i) && length(i) == 1)
#                   || (is.logical(i) && sum(i) == 1))
  singlecolumn <- ((is.integer(j) && length(j) == 1 && j > 0)
                   || (is.character(j) && length(j) == 1)
                   || (is.logical(j) && sum(j) == 1))
  if(!igiven && jgiven) {
    # x[, j] <- value
    if(singlecolumn) {
      # expecting single hypercolumn
      if(is.logical(j)) j <- names(x)[j]
      y <- get("$<-.hyperframe")(x, j, value)
    } else {
      # expecting hyperframe 
      xlist <- as.list(x)
      xlist[j] <- as.list(as.hyperframe(value))
      # the above construction accepts all indices including extra entries
      y <- do.call(hyperframe, append(xlist,
                                        list(row.names=row.names(x))))
    }
  } else {
    ## x[, ] <- value or x[i, ] <- value or x[i,j] <- value 
    ## convert indices to positive integers
    rowseq <- seq_len(dimx[1])
    colseq <- seq_len(dimx[2])
    names(rowseq) <- row.names(x)
    names(colseq) <- colnam
    I <- rowseq[i]
    J <- colseq[j]
    ## convert to lists 
    xlist <- as.list(x)
    hv <- if(is.hyperframe(value)) value else
          as.hyperframe(as.solist(value, demote=TRUE))
    vlist <- as.list(hv)
    nrowV <- dim(hv)[1]
    ncolV <- dim(hv)[2]
    if(nrowV != length(I)) {
      if(nrowV == 1) {
        ## replicate
        vlist <- lapply(vlist, rep, times=nrowV)
      } else stop(paste("Replacement value has wrong number of rows:",
                        nrowV, "should be", length(I)),
                  call.=FALSE)
    }
    if(ncolV != length(J)) {
      if(ncolV == 1) {
        ## replicate
        vlist <- rep(vlist, times=ncolV)
      } else stop(paste("Replacement value has wrong number of columns:",
                        ncolV, "should be", length(J)),
                  call.=FALSE)
    }
    ## replace entries
    for(jj in J) 
      xlist[[jj]][I] <- vlist[[jj]][I]
    ## put back together
    y <- do.call(hyperframe, append(xlist,
                                      list(row.names=row.names(x))))
  } 
  return(y)
}


split.hyperframe <- local({

  split.hyperframe <- function(x, f, drop=FALSE, ...) {
    y <- data.frame(id=seq_len(nrow(x)))
    z <- split(y, f, drop=drop)
    z <- lapply(z, getElement, name="id")
    out <- lapply(z, indexi, x=x)
    return(out)
  }

  indexi <- function(i, x) x[i,]
  
  split.hyperframe
})


"split<-.hyperframe" <- function(x, f, drop=FALSE, ..., value) {
  ix <- split(seq_len(nrow(x)), f, drop = drop, ...)
  n <- length(value)
  j <- 0
  for (i in ix) {
    j <- j%%n + 1
    x[i, ] <- value[[j]]
  }
  x
}
  
