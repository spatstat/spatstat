##
##   Math.im.R
##
##   $Revision: 1.8 $ $Date: 2020/05/08 06:09:14 $
##
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

## The following function defines what happens in Ops.im
## but the formal 'Ops' method is now in Ops.im.R

imageOp <- function(e1, e2=NULL, op) {
  ## operate on an image or pair of images
  if(is.null(e2)) {
    ## unary operation
    if(!is.element(op, c("!", "-", "+")))
      stop(paste("Unary operation", sQuote(op), "is undefined for images"),
           call.=FALSE)
    expr <- parse(text = paste(op, "e1"))
  } else {
    expr <- parse(text = paste("e1", op, "e2"))
  }
  return(do.call(eval.im, list(expr = expr)))
}                  

