#
#   distfun.R
#
#   distance function (returns a function of x,y)
#
#   $Revision: 1.20 $   $Date: 2013/12/12 05:38:59 $
#

distfun <- function(X, ...) {
  UseMethod("distfun")
}

distfun.ppp <- function(X, ..., k=1) {
  # this line forces X to be bound
  stopifnot(is.ppp(X))
  stopifnot(length(k) == 1)
  g <- function(x,y=NULL) {
    Y <- xy.coords(x, y)[c("x", "y")]
    nncross(Y, X, what="dist", k=k)
  }
  attr(g, "Xclass") <- "ppp"
  g <- funxy(g, as.rectangle(as.owin(X)))
  attr(g, "k") <- k
  class(g) <- c("distfun", class(g))
  return(g)
}

distfun.psp <- function(X, ...) {
  # this line forces X to be bound
  stopifnot(is.psp(X))
  g <- function(x,y=NULL) {
    Y <-  xy.coords(x, y)[c("x", "y")]
    nncross(Y, X, what="dist")
  }
  attr(g, "Xclass") <- "psp"
  g <- funxy(g, as.rectangle(as.owin(X)))
  class(g) <- c("distfun", class(g))
  return(g)
}

distfun.owin <- function(X, ..., invert=FALSE) {
  # this line forces X to be bound
  stopifnot(is.owin(X))
  #
  P <- as.psp(as.polygonal(X))
  #
  g <- function(x,y=NULL) {
    Y <-  xy.coords(x, y)
    inside <- inside.owin(Y$x, Y$y, X)
    D <- nncross(Y, P, what="dist")
    out <- if(!invert) ifelseAX(inside, 0, D) else ifelseXB(inside, D, 0)
    return(out)
  }
  attr(g, "Xclass") <- "owin"
  g <- funxy(g, as.rectangle(as.owin(X)))
  class(g) <- c("distfun", class(g))
  return(g)
}

as.owin.distfun <- function(W, ..., fatal=TRUE) {
  X <- get("X", envir=environment(W))
  result <- if(is.owin(X)) as.rectangle(X) else as.owin(X, ..., fatal=fatal)
  return(result)
}

as.im.distfun <- function(X, W=NULL, ...,
                           eps=NULL, dimyx=NULL, xy=NULL,
                           na.replace=NULL) {
  k <- attr(X, "k")
  if(is.null(W) && (is.null(k) || (k == 1))) {
    # use 'distmap' for speed
    env <- environment(X)
    Xdata  <- get("X",      envir=env)
    if(is.owin(Xdata)) {
      invert <- get("invert", envir=env)
      if(invert)
        Xdata <- complement.owin(Xdata)
    }
    D <- distmap(Xdata, eps=eps, dimyx=dimyx, xy=xy)
    if(!is.null(na.replace))
      D$v[is.null(D$v)] <- na.replace
  } else if(identical(attr(X, "Xclass"), "ppp")) {
    # point pattern --- use nngrid/knngrid
    env <- environment(X)
    Xdata  <- get("X",      envir=env)
    D <- nnmap(Xdata, W=W, what="dist", k=k, 
               eps=eps, dimyx=dimyx, xy=xy, na.replace=na.replace,
               ...)
  } else {
    # evaluate function at pixel centres
    D <- as.im.function(X, W=W,
                        eps=eps, dimyx=dimyx, xy=xy, na.replace=na.replace)
  }
  return(D)
}

print.distfun <- function(x, ...) {
  xtype <- attr(x, "Xclass")
  typestring <- switch(xtype,
                       ppp="point pattern",
                       psp="line segment pattern",
                       owin="window",
                       "unrecognised object")
  objname <- switch(xtype,
                    ppp="point",
                    psp="line segment",
                    "object")
  cat(paste("Distance function for", typestring, "\n"))
  X <- get("X", envir=environment(x))
  print(X)
  if(!is.null(k <- attr(x, "k")) && k > 1)
    cat(paste("Distance to", ordinal(k), "nearest", objname,
              "will be computed\n"))
  return(invisible(NULL))
}

