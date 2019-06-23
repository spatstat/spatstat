#'
#'   uniquemap.R
#'
#'   Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
#'   Licence: GNU Public Licence >= 2
#'
#'   $Revision: 1.12 $  $Date: 2019/06/23 04:41:32 $

uniquemap <- function(x) { UseMethod("uniquemap") }

uniquemap.default <- function(x) {
  result <- seqn <- seq_along(x)
  if(length(x) <= 1) return(result)
  if(is.atomic(x) && (is.factor(x) || (is.vector(x) && is.numeric(x)))) {
    if(is.factor(x)) x <- as.integer(x)
    o <- order(x, seqn)
    isfirst <- c(TRUE, (diff(x[o]) != 0))
    omap <- cumsum(isfirst)
    result <- seqn
    result[o] <- o[isfirst][omap]
    return(result)
  }
  dup <- duplicated(x)
  ux <- x[!dup]
  mapdup <- match(x[dup], ux)
  result[dup] <- which(!dup)[mapdup]
  return(result)
}

uniquemap.data.frame <- function(x) {
  n <- nrow(x)
  result <- seq_len(n)
  if(n <= 1 || !anyDuplicated(x))
    return(result)
  dup <- duplicated(x)
  uni <- which(!dup)
  for(j in which(dup)) {
    for(i in uni[uni < j]) {
      if(IdenticalRowPair(i, j, x)) {
        result[j] <- i
        break
      }
    }
  }
  return(result)
}

uniquemap.ppp <- function(x) {
  n <- npoints(x)
  seqn <- seq_len(n)
  if(n <= 1) return(seqn)
  marx <- marks(x)
  switch(markformat(marx),
         none = {
           useC <- TRUE
         },
         vector = {
           #' convert to integers if possible
           if(is.integer(marx) || is.factor(marx)) {
             marx <- as.integer(marx)
             useC <- TRUE
           } else {
             um <- unique(marx)
             if(length(um) <= 2^30) {
               marx <- match(marx, um)
               useC <- TRUE
             } else {
               useC <- FALSE
             }
           }
         },
         {
           useC <- FALSE
         })

  if(!useC) {
    #' first find duplicated spatial coordinates
    u <- uniquemap(unmark(x))
    #' add marks
    df <- cbind(data.frame(ind=seqn, uni=u), as.data.frame(marx))
    bb <- split(df, factor(u))
    #' consider each set of duplicated locations
    for(b in bb) {
      #' find unique rows of marks, as a list
      mrows <- lapply(seq_len(nrow(b)), function(i) b[i, -(1:2)])
      um <- unique(mrows)
      #' match other rows to them
      ma <- match(mrows, um)
      #' map to original index
      u[b$ind] <- b$ind[ma]
    }
    return(u)
  }

  #' unmarked or integer/factor marked
  xx <- x$x
  yy <- x$y
  o <- order(xx, seqn)

  if(is.null(marx)) {
    umap <- .C("uniqmapxy",
               n=as.integer(n),
               x=as.double(xx[o]),
               y=as.double(yy[o]),
               uniqmap=as.integer(integer(n)),
               PACKAGE="spatstat")$uniqmap
  } else {
    #' marks are (converted to) integers
    umap <- .C("uniqmap2M",
               n=as.integer(n),
               x=as.double(xx[o]),
               y=as.double(yy[o]),
               marks=as.integer(marx[o]),
               uniqmap=as.integer(integer(n)),
               PACKAGE="spatstat")$uniqmap
  }
  nodup <- (umap == 0)
  umap[nodup] <- which(nodup)
  result <- integer(n)
  result[o] <- o[umap]
  return(result)
}

uniquemap.lpp <- function(x) {
  n <- npoints(x)
  if(n <= 1 || !anyDuplicated(as.ppp(x))) return(seq_len(n))
  result <- uniquemap(as.data.frame(x))
  return(result)
}

uniquemap.ppx <- function(x) {
  uniquemap(as.data.frame(x))
}

