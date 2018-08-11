##
##   Math.im.R
##
##   $Revision: 1.7 $ $Date: 2017/01/12 03:50:22 $
##

Ops.im <- function(e1,e2=NULL){
    unary <- nargs() == 1L
    if(unary){
        if(!is.element(.Generic, c("!", "-", "+")))
            stop("Unary usage is undefined for this operation for images.")
        callstring <- paste(.Generic, "e1")
    } else {
        callstring <- paste("e1", .Generic, "e2")
    }
    expr <- parse(text = callstring)
    return(do.call(eval.im, list(expr = expr)))
}

Math.im <- function(x, ...){
    m <- do.call(.Generic, list(x$v, ...))
    rslt <- im(m, xcol = x$xcol, yrow = x$yrow, xrange = x$xrange,
               yrange = x$yrange, unitname = unitname(x))
    return(rslt)
}

Summary.im <- function(..., na.rm=FALSE, drop=TRUE){
  argh <- list(...)
  ims <- sapply(argh, is.im)
  argh[ims] <- lapply(argh[ims], getElement, name="v")
  do.call(.Generic, c(argh, list(na.rm = na.rm || drop)))
}

Complex.im <- function(z){
    m <- do.call(.Generic, list(z=z$v))
    rslt <- im(m, xcol = z$xcol, yrow = z$yrow, xrange = z$xrange,
               yrange = z$yrange, unitname = unitname(z))
    return(rslt)
}
