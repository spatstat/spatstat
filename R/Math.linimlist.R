##
##   Math.linimlist.R
##
##   $Revision: 1.1 $ $Date: 2020/10/31 08:47:12 $
##

Math.linimlist <- function(x, ...){
  solapply(x, .Generic, ...)
}

Complex.linimlist <- function(z){
  solapply(z, .Generic)
}

Summary.linimlist <- function(..., na.rm=TRUE){
  argh <- expandSpecialLists(list(...))
  if(length(names(argh)) > 0) {
    isim <- sapply(argh, is.im)
    names(argh)[isim] <- ""
  }
  do.call(.Generic, c(argh, list(na.rm=na.rm)))
}

#'   Due to the dispatch mechanism, Ops.linim and Ops.linimlist must be identical
#'   if we want to handle combinations of linimlist and im.
#'   (See 'Math.linim.R' for the definition of 'LinimOp')

Ops.linimlist <- Ops.linim <- function(e1,e2=NULL){ LinimListOp(e1, e2, .Generic) }

LinimListOp <- function(e1, e2=NULL, op) {
  if(is.null(e2)) {
    #' unary operation
    result <- if(is.im(e1)) LinimOp(e1, op=op) else solapply(e1, LinimOp, op=op)
    return(result)
  } 
  #' binary operation
  single1 <- !inherits(e1, c("linimlist", "solist"))
  single2 <- !inherits(e2, c("linimlist", "solist"))
  if(single1 && single2) return(LinimOp(e1, e2, op))
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
  v <- mapply(LinimOp, e1=unname(e1list), e2=unname(e2list),
              MoreArgs=list(op=op),
              SIMPLIFY=FALSE)
  names(v) <- outnames
  return(as.solist(v))
}

