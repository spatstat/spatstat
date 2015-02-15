##
##   Math.linim.R
##
##   $Revision: 1.3 $ $Date: 2015/02/15 10:50:01 $
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
    rslt <- linim(L, Z, df=df)
    return(rslt)
}

Summary.linim <- function(..., na.rm){
    args <- list(...)
    argp <- lapply(args, "[")
    argd <- if(is.element(.Generic, c("sum", "prod"))) list() else 
            lapply(lapply(args, attr, which="df"), getElement, name="values")
    do.call(.Generic, c(argp, argd, na.rm = na.rm))
}

Complex.linim <- function(z){
    L <- attr(z, "L")
    df <- attr(z, "df")
    m <- do.call(.Generic, list(z=z[drop=TRUE]))
    Z <- im(m, xcol = z$xcol, yrow = z$yrow, xrange = z$xrange,
               yrange = z$yrange, unitname = unitname(z))
    df$values <- do.call(.Generic, list(z=df$values))
    rslt <- linim(L, Z, df=df)
    return(rslt)
}
