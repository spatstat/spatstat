#
#   nnfun.R
#
#   nearest neighbour function (returns a function of x,y)
#
#   $Revision: 1.8 $   $Date: 2019/09/16 10:14:28 $
#

nnfun <- function(X, ...) {
  UseMethod("nnfun")
}

nnfun.ppp <- function(X, ..., k=1, value=c("index", "mark")) {
  # this line forces X to be bound
  stopifnot(is.ppp(X))
  if(length(k) != 1) stop("k should be a single integer")
  value <- match.arg(value)
  switch(value,
         index = {
           gi <- function(x,y=NULL) {
             Y <- xy.coords(x, y)[c("x", "y")]
             nncross(Y, X, what="which", k=k)
           }
           attr(gi, "Xclass") <- "ppp"
           g <- funxy(gi, as.rectangle(as.owin(X)))
         },
         mark = {
           stopifnot(is.marked(X))
           marx <- as.data.frame(marks(X))[,1]
           gm <- function(x,y=NULL) {
             Y <- xy.coords(x, y)[c("x", "y")]
             marx[nncross(Y, X, what="which", k=k)]
           }
           attr(gm, "Xclass") <- "ppp"
           g <- funxy(gm, as.rectangle(as.owin(X)))
         })
  class(g) <- c("nnfun", class(g))
  return(g)
}

nnfun.psp <- function(X, ..., value=c("index", "mark")) {
  # this line forces X to be bound
  stopifnot(is.psp(X))
  value <- match.arg(value)
  switch(value,
         index = {
           gi <- function(x,y=NULL) {
             Y <-  xy.coords(x, y)[c("x", "y")]
             nncross(Y, X, what="which")
           }
           attr(gi, "Xclass") <- "psp"
           g <- funxy(gi, as.rectangle(as.owin(X)))
         },
         mark = {
           stopifnot(is.marked(X))
           marx <- as.data.frame(marks(X))[,1]
           gm <- function(x,y=NULL) {
             Y <-  xy.coords(x, y)[c("x", "y")]
             marx[nncross(Y, X, what="which")]
           }
           attr(gm, "Xclass") <- "psp"
           g <- funxy(gm, as.rectangle(as.owin(X)))
         })
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
                           na.replace=NULL, approx=TRUE) {
  if(approx && is.null(W)) {
    env <- environment(X)
    Xdata  <- get("X", envir=env)
    k <- mget("k", envir=env, inherits=FALSE, ifnotfound=list(1))[[1L]]
    value <- mget("value", envir=env, ifnotfound=list(NULL))[[1L]] 
    Z <- nnmap(Xdata, k=k, what="which", eps=eps, dimyx=dimyx, xy=xy)
    if(identical(value, "mark")) {
      marx <- get("marx", envir=env)
      Z <- eval.im(marx[Z])
    }
    if(!is.null(na.replace))
      Z$v[is.null(Z$v)] <- na.replace
    return(Z)
  }
  if(is.null(W)) W <- Window(X)
  result <- as.im.function(X, W=W,
                           eps=eps, dimyx=dimyx, xy=xy,
                           na.replace=na.replace, ...)
  return(result)
}

print.nnfun <- function(x, ...) {
  env <- environment(x)
  X <- get("X", envir=env)
  k <- mget("k", envir=env, inherits=FALSE, ifnotfound=list(1))[[1L]]
  v <- mget("value", envir=env, ifnotfound=list(NULL))[[1L]]
  xtype <- attr(x, "Xclass")
  typestring <- switch(xtype,
                       ppp="point pattern",
                       psp="line segment pattern",
                       paste("object of class", sQuote(xtype)))
  Kth <- if(k == 1) "Nearest" else paste0(ordinal(k), "-Nearest")
  cat(paste(Kth, "Neighbour",
            if(is.null(v)) "Index" else "Mark",
            "function for ", typestring, "\n"))
  print(X)

  return(invisible(NULL))
}

