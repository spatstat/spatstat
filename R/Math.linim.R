##
##   Math.linim.R
##
##   $Revision: 1.8 $ $Date: 2020/04/09 01:59:33 $
##

Ops.linim <- function(e1,e2=NULL){
    unary <- nargs() == 1L
    if(unary){
        if(!is.element(.Generic, c("!", "-", "+")))
            stop("Unary usage is undefined for this operation for images.")
        callstring <- paste(.Generic, "e1")
    } else {
        callstring <- paste("e1", .Generic, "e2")
    }
    expr <- parse(text = callstring)
    return(do.call(eval.linim, list(expr = expr)))
}

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
