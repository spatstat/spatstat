#
#   funxy.R
#
#   Class of functions of x,y location with a spatial domain
#
#   $Revision: 1.2 $   $Date: 2012/10/13 09:35:28 $
#

funxy <- function(f, W=NULL) {
  stopifnot(is.function(f))
  stopifnot(is.owin(W))
  class(f) <- c("funxy", class(f))
  attr(f, "W") <- W
  return(f)
}

print.funxy <- function(x, ...) {
  cat(paste("function(x,y) of class", sQuote("funxy"), "\n"))
  print(as.owin(x))
}

as.owin.funxy <- function(W, ..., fatal=TRUE) {
  W <- attr(W, "W")
  as.owin(W, ..., fatal=fatal)
}

domain.funxy <- Window.funxy <- function(X, ...) { as.owin(X) }

#   Note that 'distfun' (and other classes inheriting from funxy)
#   has a method for as.owin that takes precedence over as.owin.funxy
#   and this will affect the behaviour of the following plot methods
#   because 'distfun' does not have its own plot method.

plot.funxy <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  W <- as.owin(x)
  do.call("do.as.im",
          resolve.defaults(list(x, action="plot"),
                           list(...),
                           list(main=xname, W=W)))
  invisible(NULL)
}

contour.funxy <- function(x, ...) {
  xname <- deparse(substitute(x))
  W <- as.owin(x)
  do.call("do.as.im",
          resolve.defaults(list(x, action="contour"),
                           list(...),
                           list(main=xname, W=W)))
  invisible(NULL)
}

persp.funxy <- function(x, ...) {
  xname <- deparse(substitute(x))
  W <- as.rectangle(as.owin(x))
  do.call("do.as.im",
          resolve.defaults(list(x, action="persp"),
                           list(...),
                           list(main=xname, W=W)))
  invisible(NULL)
}

