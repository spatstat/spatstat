##
##   Math.imlist.R
##
##   $Revision: 1.2 $ $Date: 2017/01/07 02:40:15 $
##

Math.imlist <- function(x, ...){
  solapply(x, .Generic, ...)
}

Complex.imlist <- function(z){
  solapply(z, .Generic)
}

Summary.imlist <- function(..., na.rm=TRUE){
  argh <- expandSpecialLists(list(...))
  if(length(names(argh)) > 0) {
    isim <- sapply(argh, is.im)
    names(argh)[isim] <- ""
  }
  do.call(.Generic, c(argh, list(na.rm=na.rm)))
}

Ops.imlist <- function(e1,e2=NULL){
  if(nargs() == 1L) {
    #' unary operation
    solapply(e1, .Generic)
  } else {
    #' binary operation
    if(inherits(e2, "imlist")) {
      as.solist(mapply(.Generic, unname(e1), unname(e2), SIMPLIFY=FALSE))
    } else {
      solapply(e1, .Generic, e2=e2)
    }
  }
}


