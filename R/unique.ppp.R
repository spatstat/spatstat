#
#   unique.ppp.R
#
# $Revision: 1.33 $  $Date: 2018/11/02 01:19:08 $
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
  if(rule == "unmark")
    x <- unmark(x)
  n <- npoints(x)
  switch(markformat(x),
         none = {
           # unmarked points
           # check for duplication of x and y separately (a necessary condition)
           xx <- x$x
           yy <- x$y
           possible <- duplicated(xx) & duplicated(yy)
           if(!any(possible))
             return(possible)
           # split by x coordinate of duplicated x values
           result <- possible
           xvals <- unique(xx[possible])
           for(xvalue in xvals) {
             sub <- (xx == xvalue)
           # compare y values
             result[sub] <- duplicated(yy[sub])
           }
         },
         vector = {
           # marked points - split by mark value
           m <- marks(x)
           um <- if(is.factor(m)) levels(m) else unique(m)
           xx <- unmark(x)
           result <- logical(n)
           for(i in seq_along(um)) {
             sub <- (m == um[i])
             result[sub] <- duplicated.ppp(xx[sub])
           }
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
  anyDuplicated(as.data.frame(x), ...)
}

## utility to check whether two rows are identical

IdenticalRows <- local({
  id <- function(i,j, a, b=a) {
    ai <- a[i,]
    bj <- b[j,]
    row.names(ai) <- row.names(bj) <- NULL
    identical(ai, bj)
  }
  Vectorize(id, c("i", "j"))
})
    

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


