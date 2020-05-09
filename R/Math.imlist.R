##
##   Math.imlist.R
##
##   $Revision: 1.6 $ $Date: 2020/05/09 03:32:49 $
##
##  'Ops' method is now in Ops.im.R

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

#'   Due to the dispatch mechanism, Ops.im and Ops.imlist must be identical
#'   if we want to handle combinations of imlist and im.
#'   (See 'Math.im.R' for the definition of 'imageOp')

Ops.imlist <- Ops.im <- function(e1,e2=NULL){ imagelistOp(e1, e2, .Generic) }

imagelistOp <- function(e1, e2=NULL, op) {
  if(is.null(e2)) {
    #' unary operation
    result <- if(is.im(e1)) imageOp(e1, op=op) else solapply(e1, imageOp, op=op)
    return(result)
  } 
  #' binary operation
  single1 <- !inherits(e1, c("imlist", "solist"))
  single2 <- !inherits(e2, c("imlist", "solist"))
  if(single1 && single2) return(imageOp(e1, e2, op))
  if(single1 && !single2) {
    e1list <- rep(list(e1), length(e2))
    e2list <- e2
    outnames <- names(e2)
  } else if(!single1 && single2) {
    e1list <- e1
    e2list <- rep(list(e2), length(e1))
    outnames <- names(e1)
  } else {
    e1list <- e1
    e2list <- e2
    if(length(e1) != length(e2))
      stop(paste("Lists of images have incompatible lengths:",
                 length(e1), "!=", length(e2)),
           call.=FALSE)
    outnames <- names(e1) %orifnull% names(e2)
  }
  #' compute
  v <- mapply(imageOp, e1=unname(e1list), e2=unname(e2list),
              MoreArgs=list(op=op),
              SIMPLIFY=FALSE)
  names(v) <- outnames
  return(as.solist(v))
}


