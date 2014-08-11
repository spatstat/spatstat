#
#   distfun.R
#
#   distance function (returns a function of x,y)
#
#   $Revision: 1.15 $   $Date: 2013/04/25 06:37:43 $
#

distfun <- function(X, ...) {
  UseMethod("distfun")
}

distfun.ppp <- function(X, ...) {
  # this line forces X to be bound
  stopifnot(is.ppp(X))
  g <- function(x,y=NULL) {
    Y <- xy.coords(x, y)[c("x", "y")]
    nncross(Y, X, what="dist")
  }
  attr(g, "Xclass") <- "ppp"
  g <- funxy(g, as.rectangle(as.owin(X)))
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
  if(X$type == "mask" && (!(spatstat.options("gpclib") && require(gpclib)))) {
    warning("Polygon calculations unavailable; using distmap")
    discrete <- TRUE
    D <- if(!invert) distmap(X) else distmap(complement.owin(X))
  } else {
    discrete <- FALSE
    P <- as.psp(as.polygonal(X))
  }
  g <- function(x,y=NULL) {
    Y <-  xy.coords(x, y)
    if(discrete)
      return(D[Y])
    inside <- inside.owin(Y$x, Y$y, X)
    D <- nncross(Y, P, what="dist")
    out <- if(!invert) ifelse(inside, 0, D) else ifelse(inside, D, 0)
    return(out)
  }
  attr(g, "Xclass") <- "owin"
  g <- funxy(g, as.rectangle(as.owin(X)))
  class(g) <- c("distfun", class(g))
  return(g)
}

as.owin.distfun <- function(W, ..., fatal=TRUE) {
  X <- get("X", envir=environment(W))
  as.owin(X, ..., fatal=fatal)
}

as.im.distfun <- function(X, W=NULL, ...,
                           eps=NULL, dimyx=NULL, xy=NULL,
                           na.replace=NULL) {
  if(is.null(W)) {
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
    return(D)
  }
  # use as.im.function
  NextMethod("as.im")
}

print.distfun <- function(x, ...) {
  xtype <- attr(x, "Xclass")
  typestring <- switch(xtype,
                       ppp="point pattern",
                       psp="line segment pattern",
                       owin="window",
                       "unrecognised object")
  cat(paste("Distance function for", typestring, "\n"))
  X <- get("X", envir=environment(x))
  print(X)
  return(invisible(NULL))
}

