#'
#'   uniquemap.R
#'
#'   Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
#'   Licence: GNU Public Licence >= 2
#'
#'   $Revision: 1.10 $  $Date: 2019/05/21 09:23:35 $

uniquemap <- function(x) { UseMethod("uniquemap") }

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

