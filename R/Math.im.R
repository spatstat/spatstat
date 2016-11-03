##
##   Math.im.R
##
##   $Revision: 1.4 $ $Date: 2016/11/03 09:11:24 $
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
    m <- do.call(.Generic, list(x[,,drop=FALSE], ...))
    rslt <- im(m, xcol = x$xcol, yrow = x$yrow, xrange = x$xrange,
               yrange = x$yrange, unitname = unitname(x))
    return(rslt)
}

Summary.im <- function(..., na.rm){
  argh <- list(...)
  ims <- sapply(argh, is.im)
  argh[ims] <- lapply(argh[ims], "[")
  do.call(.Generic, argh)
}

Complex.im <- function(z){
    m <- do.call(.Generic, list(z=z[drop=TRUE]))
    rslt <- im(m, xcol = z$xcol, yrow = z$yrow, xrange = z$xrange,
               yrange = z$yrange, unitname = unitname(z))
    return(rslt)
}
