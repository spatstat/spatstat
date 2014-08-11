#
#   nnfun.R
#
#   nearest neighbour function (returns a function of x,y)
#
#   $Revision: 1.1 $   $Date: 2012/10/13 10:09:55 $
#

nnfun <- function(X, ...) {
  UseMethod("nnfun")
}

nnfun.ppp <- function(X, ...) {
  # this line forces X to be bound
  stopifnot(is.ppp(X))
  g <- function(x,y=NULL) {
    Y <- xy.coords(x, y)[c("x", "y")]
    nncross(Y, X, what="which")
  }
  attr(g, "Xclass") <- "ppp"
  g <- funxy(g, as.rectangle(as.owin(X)))
  class(g) <- c("nnfun", class(g))
  return(g)
}

nnfun.psp <- function(X, ...) {
  # this line forces X to be bound
  stopifnot(is.psp(X))
  g <- function(x,y=NULL) {
    Y <-  xy.coords(x, y)[c("x", "y")]
    nncross(Y, X, what="which")
  }
  attr(g, "Xclass") <- "psp"
  g <- funxy(g, as.rectangle(as.owin(X)))
  class(g) <- c("nnfun", class(g))
  return(g)
}

as.owin.nnfun <- function(W, ..., fatal=TRUE) {
  X <- get("X", envir=environment(W))
  as.owin(X, ..., fatal=fatal)
}

as.im.nnfun <- function(X, W=NULL, ...,
                           eps=NULL, dimyx=NULL, xy=NULL,
                           na.replace=NULL) {
  if(is.null(W)) {
    # use 'distmap' for speed
    env <- environment(X)
    Xdata  <- get("X",      envir=env)
    D <- distmap(Xdata, eps=eps, dimyx=dimyx, xy=xy)
    Z <- attr(D, "index") 
    if(!is.null(na.replace))
      Z$v[is.null(Z$v)] <- na.replace
    return(Z)
  }
  # use as.im.function
  NextMethod("as.im")
}

print.nnfun <- function(x, ...) {
  xtype <- attr(x, "Xclass")
  typestring <- switch(xtype,
                       ppp="point pattern",
                       psp="line segment pattern",
                       owin="window",
                       "unrecognised object")
  cat(paste("Nearest Neighbour Index function for", typestring, "\n"))
  X <- get("X", envir=environment(x))
  print(X)
  return(invisible(NULL))
}

