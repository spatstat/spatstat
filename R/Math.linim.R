##
##   Math.linim.R
##
##   $Revision: 1.9 $ $Date: 2020/10/31 08:47:06 $
##

Math.linim <- function(x, ...){
    m <- do.call(.Generic, list(x[,,drop=FALSE], ...))
    Z <- im(m, xcol = x$xcol, yrow = x$yrow, xrange = x$xrange,
            yrange = x$yrange, unitname = unitname(x))
    df <- attr(x, "df")
    df$values <- do.call(.Generic, list(df$values, ...))
    L <- attr(x, "L")
    rslt <- linim(L, Z, df=df, restrict=FALSE)
    return(rslt)
}

Summary.linim <- function(..., na.rm, finite){
  if(missing(finite)) finite <- FALSE
  if(missing(na.rm)) na.rm <- FALSE
  argh <- list(...)
  values <- lapply(argh, "[")
  dfvalues <- if(is.element(.Generic, c("sum", "prod"))) list() else 
              lapply(lapply(argh, attr, which="df"), getElement, name="values")
  vals <- unlist(c(values, dfvalues))
  logique <- is.element(.Generic, c("all", "any"))
  vals <- if(logique) as.logical(vals) else as.numeric(vals)
  if(finite && !logique) {
    vals <- vals[is.finite(vals)]
  } else if(na.rm) {
    vals <- vals[!is.na(vals)]
  }
  do.call(.Generic, list(vals))
}

Complex.linim <- function(z){
    L <- attr(z, "L")
    df <- attr(z, "df")
    m <- do.call(.Generic, list(z=z[,,drop=FALSE]))
    Z <- im(m, xcol = z$xcol, yrow = z$yrow, xrange = z$xrange,
               yrange = z$yrange, unitname = unitname(z))
    df$values <- do.call(.Generic, list(z=df$values))
    rslt <- linim(L, Z, df=df, restrict=FALSE)
    return(rslt)
}

## This function defines what happens inside 'Ops.linim'
## The formal definition of 'Ops.linim' is now in 'Math.linimlist.R'

LinimOp <- function(e1, e2=NULL, op) {
  ## operate on a linim or pair of linim with fallback to im
  if(is.null(e2)) {
    ## unary operation
    if(!is.element(op, c("!", "-", "+")))
      stop(paste("Unary operation", sQuote(op), "is undefined for images"),
           call.=FALSE)
    expr <- parse(text = paste(op, "e1"))
    netted <- is.linim(e1)
  } else {
    expr <- parse(text = paste("e1", op, "e2"))
    netted <- is.linim(e1) && is.linim(e2)
  }
  result <- if(netted) {
              do.call(eval.linim, list(expr=expr))
            } else {
              do.call(eval.im, list(expr = expr))
            }
  return(result)
}                  

