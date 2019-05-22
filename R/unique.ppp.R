#
#   unique.ppp.R
#
# $Revision: 1.37 $  $Date: 2019/05/22 09:04:57 $
#
# Methods for 'multiplicity' co-authored by Sebastian Meyer
# Copyright 2013 Adrian Baddeley and Sebastian Meyer 

unique.ppp <- function(x, ..., warn=FALSE) {
  verifyclass(x, "ppp")
  dupe <- duplicated.ppp(x, ...) 
  if(!any(dupe)) return(x)
  if(warn) warning(paste(sum(dupe), "duplicated points were removed"),
                   call.=FALSE)
  return(x[!dupe])
}

duplicated.ppp <- function(x, ...,
                           rule=c("spatstat", "deldir", "unmark")) {
  verifyclass(x, "ppp")
  rule <- match.arg(rule)
  if(rule == "deldir")
    return(deldir::duplicatedxy(x))
  n <- npoints(x)
  xloc <- unmark(x)
  if(!anyDuplicated(xloc))
    return(logical(n)) # i.e. vector of FALSE
  if(rule == "unmark")
    x <- xloc
  switch(markformat(x),
         none = {
           #' unmarked points
           u <- uniquemap(x)
           result <- (u != seq_along(u))
         },
         vector = {
           #' marked points - convert mark to integer
           m <- marks(x)
           if(is.factor(m)) {
             marks(x) <- as.integer(m)
           } else {
             um <- unique(m)
             marks(x) <- match(m, um)
           }
           result <- duplicated(as.data.frame(x))
         },
         dataframe = {
           result <- duplicated(as.data.frame(x))
         },
         # the following are currently not supported
         hyperframe = {
           result <- duplicated(as.data.frame(x))
         }, 
         list = {
           result <- duplicated(as.data.frame(as.hyperframe(x)))
         },
         stop(paste("Unknown mark type", sQuote(markformat(x))))
         )
  return(result)
}

anyDuplicated.ppp <- function(x, ...) {
  #' first check duplication of coordinates using fast code
  n <- npoints(x)
  if(n <= 1) return(FALSE)
  xx <- x$x
  yy <- x$y
  o <- order(xx, seq_len(n))
  anydupXY <- .C("anydupxy",
                 n=as.integer(n),
                 x=as.double(xx[o]),
                 y=as.double(yy[o]),
                 anydup=as.integer(integer(1)),
                 PACKAGE="spatstat")$anydup
  anydupXY && (!is.marked(x) || anyDuplicated(as.data.frame(x), ...))
}

## utility to check whether two rows are identical

IdenticalRowPair <- function(i,j, a, b=a) {
  #' i and j are row indices (single integers)
  ai <- a[i,]
  bj <- b[j,]
  row.names(ai) <- row.names(bj) <- NULL
  identical(ai, bj)
}

## vectorised

IdenticalRows <- function(i, j, a, b=a) {
  #' i and j are row index vectors of equal length
  #' result[k] = identical( a[i[k],]  , b[j[k],] )
  Mo <- if(missing(b)) list(a=a) else list(a=a, b=b)
  mapply(IdenticalRowPair, i=i, j=j, MoreArgs=Mo,
         SIMPLIFY=TRUE, USE.NAMES=FALSE)
}

## .......... multiplicity .............

multiplicity <- function(x) {
  UseMethod("multiplicity")
}
  
multiplicity.ppp <- function(x) {
  verifyclass(x, "ppp")
  np <- npoints(x)
  if(np == 0) return(integer(0))
  cl <- closepairs(x, 0, what="indices")
  I <- cl$i
  J <- cl$j
  if(length(I) == 0)
    return(rep.int(1L, np))
  switch(markformat(x),
         none = { },
         vector = {
           marx <- as.data.frame(marks(x))
           agree <- IdenticalRows(I, J, marx)
           I <- I[agree]
           J <- J[agree]
         },
         dataframe = {
           marx <- marks(x)
           agree <- IdenticalRows(I, J, marx)
           I <- I[agree]
           J <- J[agree]
         },
         hyperframe = {
           marx <- as.data.frame(marks(x)) # possibly discards columns
           agree <- IdenticalRows(I, J, marx)
           I <- I[agree]
           J <- J[agree]
         }, 
         list = stop("Not implemented for lists of marks")
         )
  if(length(I) == 0)
    return(rep.int(1L, np))
  JbyI <- split(J, factor(I, levels=1:np))
  result <- 1 + lengths(JbyI)
  return(result)
}
  
multiplicity.data.frame <- function (x) {
  if(all(unlist(lapply(x, is.numeric))))
    return(multiplicityNumeric(as.matrix(x)))
  ## result template (vector of 1's)
  result <- setNames(rep.int(1L, nrow(x)), rownames(x))
  ## check for duplicates (works for data frames, arrays and vectors)
  ## CAVE: comparisons are based on a character representation of x
  if (!any(dup <- duplicated(x)))
    return(result)
  ux <- x[!dup, , drop=FALSE]
  dx <- x[dup,  , drop=FALSE]
  nu <- nrow(ux)
  nd <- nrow(dx)
  hit <- outer(seq_len(nu), seq_len(nd), IdenticalRows, a=ux, b=dx)
  counts <- as.integer(1L + .rowSums(hit, nu, nd))
  result[!dup] <- counts
  dumap <- apply(hit, 2, match, x=TRUE) # equivalent to min(which(z))
  result[dup] <- counts[dumap]
  return(result)
}

### multiplicity method for NUMERIC arrays, data frames, and vectors
### This implementation is simply based on checking for dist(x)==0

multiplicityNumeric <- function(x)
{
  if (anyDuplicated(x)) {
    distmat <- as.matrix(dist(x, method="manhattan"))  # faster than euclid.
    result <- as.integer(rowSums(distmat == 0))        # labels are kept
    if(is.null(names(result)))
      names(result) <- seq_along(result)
  } else {                                             # -> vector of 1's
    nx <- NROW(x)
    labels <- if (length(dim(x))) rownames(x) else names(x)
    if (is.null(labels)) labels <- seq_len(nx)
    result <- setNames(rep.int(1L, nx), labels)
  }
  return(result)
}

### multiplicity method for arrays, data frames, and vectors (including lists)
### It also works for non-numeric data, since it is based on duplicated().

multiplicity.default <- function (x) {
  if(is.numeric(x))
    return(multiplicityNumeric(x))
  nx <- NROW(x)                   # also works for a vector x
  ## result template (vector of 1's)
  labels <- if (length(dim(x))) rownames(x) else names(x)
  if (is.null(labels)) labels <- seq_len(nx)
  result <- setNames(rep.int(1L, nx), labels)
  ## check for duplicates (works for data frames, arrays and vectors)
  ## CAVE: comparisons are based on a character representation of x
  if (!any(dup <- duplicated(x)))
    return(result)

  ## convert x to a matrix for IdenticalRows()
  x <- as.matrix(x)
  dimnames(x) <- NULL             # discard any names!
  ux <- x[!dup, , drop=FALSE]
  dx <- x[dup,  , drop=FALSE]
  nu <- nrow(ux)
  nd <- nrow(dx)
  hit <- outer(seq_len(nu), seq_len(nd), IdenticalRows, a=ux, b=dx)
  counts <- as.integer(1L + .rowSums(hit, nu, nd))
  dumap <- apply(hit, 2, match, x=TRUE) # was: function(z) min(which(z)))
  result[!dup] <- counts
  result[dup]  <- counts[dumap]
  return(result)
}


