#
#   nnfun.R
#
#   nearest neighbour function (returns a function of x,y)
#
#   $Revision: 1.4 $   $Date: 2013/12/11 09:22:56 $
#

nnfun <- function(X, ...) {
  UseMethod("nnfun")
}

nnfun.ppp <- function(X, ..., k=1) {
  # this line forces X to be bound
  stopifnot(is.ppp(X))
  if(length(k) != 1) stop("k should be a single integer")
  g <- function(x,y=NULL) {
    Y <- xy.coords(x, y)[c("x", "y")]
    nncross(Y, X, what="which", k=k)
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

domain.nnfun <- Window.nnfun <- function(X, ...) { as.owin(X) }

as.im.nnfun <- function(X, W=NULL, ...,
                           eps=NULL, dimyx=NULL, xy=NULL,
                           na.replace=NULL) {
  if(is.null(W)) {
    env <- environment(X)
    Xdata  <- get("X", envir=env)
    k <- mget("k", envir=env, inherits=FALSE, ifnotfound=list(1))[[1]]
    Z <- nnmap(Xdata, k=k, what="which", eps=eps, dimyx=dimyx, xy=xy)
    if(!is.null(na.replace))
      Z$v[is.null(Z$v)] <- na.replace
    return(Z)
  }
  # use as.im.function
  NextMethod("as.im")
}

print.nnfun <- function(x, ...) {
  env <- environment(x)
  X <- get("X", envir=env)
  k <- mget("k", envir=env, inherits=FALSE, ifnotfound=list(1))[[1]]
  xtype <- attr(x, "Xclass")
  typestring <- switch(xtype,
                       ppp="point pattern",
                       psp="line segment pattern",
                       paste("object of class", sQuote(xtype)))
  Kth <- if(k == 1) "Nearest" else paste0(ordinal(k), "-Nearest")
  cat(paste(Kth, "Neighbour Index function for ", typestring, "\n"))
  print(X)
  return(invisible(NULL))
}

