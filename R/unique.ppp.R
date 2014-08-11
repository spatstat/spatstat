#
#   unique.ppp.R
#
# $Revision: 1.21 $  $Date: 2013/11/19 03:16:20 $
#

unique.ppp <- function(x, ..., warn=FALSE) {
  verifyclass(x, "ppp")
  dupe <- duplicated.ppp(x, ...)
  if(!any(dupe)) return(x)
  if(warn) warning(paste(sum(dupe), "duplicated points were removed"),
                   call.=FALSE)
  return(x[!dupe])
}

duplicated.ppp <- function(x, ..., rule=c("spatstat", "deldir")) {
  verifyclass(x, "ppp")
  rule <- match.arg(rule)
  if(rule == "deldir")
    return(duplicatedxy(x))
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
         listof = {
           result <- duplicated(as.data.frame(as.hyperframe(x)))
         },
         stop(paste("Unknown mark type", sQuote(markformat(x))))
         )
  return(result)
}

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
    return(rep.int(1, np))
  switch(markformat(x),
         none = { },
         vector = {
           marx <- marks(x)
           agree <- (marx[I] == marx[J])
           I <- I[agree]
           J <- J[agree]
         },
         dataframe = {
           marx <- marks(x)
           agree <- apply(marx[I, ,drop=FALSE] == marx[J, ,drop=FALSE], 1, all)
           I <- I[agree]
           J <- J[agree]
         },
         hyperframe =, 
         listof = stop("Not implemented for hyperframes or lists of marks")
         )
  if(length(I) == 0)
    return(rep.int(1, np))
  JbyI <- split(J, factor(I, levels=1:np))
  result <- 1 + sapply(JbyI, length)
  return(result)
}
  
multiplicity.data.frame <- local({

  id <- function(i,j, a, b) identical(a[i,], b[j,])
  IdenticalRows <- Vectorize(id, c("i", "j")) 

  multiplicity.data.frame <- function(x) {
    dup <- duplicated(x)
    nx <- nrow(x)
    if(!any(dup))
      return(rep.int(1, nx))
    ux <- x[!dup, , drop=FALSE]
    dx <- x[dup,  , drop=FALSE]
    row.names(ux) <- NULL
    row.names(dx) <- NULL
    nu <- nrow(ux)
    nd <- nrow(dx)
    hit <- outer(seq_len(nu), seq_len(nd), IdenticalRows, a=ux, b=dx)
    counts <- 1 + rowSums(hit)
    result <- numeric(nx)
    result[!dup] <- counts
    dumap <- apply(hit, 2, function(z) min(which(z)))
    result[dup] <- counts[dumap]
    return(result)
  }

  multiplicity.data.frame
})

