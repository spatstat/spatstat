#'
#'   uniquemap.R
#'
#'   Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2019
#'   Licence: GNU Public Licence >= 2
#'
#'   $Revision: 1.5 $  $Date: 2019/05/17 09:14:23 $

uniquemap <- function(x) { UseMethod("uniquemap") }

uniquemap.ppp <- function(x) {
  n <- npoints(x)
  seqn <- seq_len(n)
  if(n <= 1) return(seqn)
  xx <- x$x
  yy <- x$y
  o <- order(xx, seqn)
  umap <- .C("uniqmapxy",
          n=as.integer(n),
          x=as.double(xx[o]),
          y=as.double(yy[o]),
          uniqmap=as.integer(integer(n)),
          PACKAGE="spatstat")$uniqmap
  nodup <- (umap == 0)
  umap[nodup] <- which(nodup)
  result <- integer(n)
  result[o] <- o[umap]
  return(result)
}




