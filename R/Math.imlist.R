##
##   Math.imlist.R
##
##   $Revision: 1.4 $ $Date: 2017/08/15 03:46:57 $
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
    return(solapply(e1, .Generic))
  } 
  #' binary operation
  if(inherits(e2, "imlist")) {
    #' two image lists - must have equal length
    v <- mapply(.Generic, unname(e1), unname(e2), SIMPLIFY=FALSE)
    names(v) <- names(e1)
    return(as.solist(v))
  }
  #' other binary operation e.g. imlist + constant, imlist + im
  return(solapply(e1, .Generic, e2=e2))
}


