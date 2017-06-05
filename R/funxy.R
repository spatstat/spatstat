#
#   funxy.R
#
#   Class of functions of x,y location with a spatial domain
#
#   $Revision: 1.14 $   $Date: 2017/06/05 10:31:58 $
#

spatstat.xy.coords <- function(x,y) {
  if(missing(y) || is.null(y)) {
    xy <- if(is.ppp(x) || is.lpp(x)) coords(x) else
          if(checkfields(x, c("x", "y"))) x else 
          stop("Argument y is missing", call.=FALSE)
    x <- xy$x
    y <- xy$y
  }
  xy.coords(x,y)[c("x","y")]
}

funxy <- function(f, W=NULL) {
  stopifnot(is.function(f))
  stopifnot(is.owin(W))
  if(!identical(names(formals(f))[1:2], c("x", "y")))
    stop("The first two arguments of f should be named x and y", call.=FALSE)
  if(is.primitive(f))
    stop("Not implemented for primitive functions", call.=FALSE)
  ## copy 'f' including formals, environment, attributes
  h <- f
  ## make new function body:
  ## paste body of 'f' into last line of 'spatstat.xy.coords'
  b <- body(spatstat.xy.coords)
  b[[length(b)]] <- body(f)
  ## transplant the body 
  body(h) <- b
  ## reinstate attributes
  attributes(h) <- attributes(f)
  ## stamp it
  class(h) <- c("funxy", class(h))
  attr(h, "W") <- W
  attr(h, "f") <- f
  return(h)  
}

print.funxy <- function(x, ...) {
  nama <- names(formals(x))
  splat(paste0("function", paren(paste(nama,collapse=","))),
        "of class", sQuote("funxy"))
  print(as.owin(x))
  splat("\nOriginal function definition:")
  print(attr(x, "f"))
}

summary.funxy <- function(object, ...) { print(object, ...) }

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
  do.call(do.as.im,
          resolve.defaults(list(x, action="plot"),
                           list(...),
                           list(main=xname, W=W)))
  invisible(NULL)
}

contour.funxy <- function(x, ...) {
  xname <- deparse(substitute(x))
  W <- as.owin(x)
  do.call(do.as.im,
          resolve.defaults(list(x, action="contour"),
                           list(...),
                           list(main=xname, W=W)))
  invisible(NULL)
}

persp.funxy <- function(x, ...) {
  xname <- deparse(substitute(x))
  W <- as.rectangle(as.owin(x))
  do.call(do.as.im,
          resolve.defaults(list(x, action="persp"),
                           list(...),
                           list(main=xname, W=W)))
  invisible(NULL)
}

hist.funxy <- function(x, ..., xname) {
  if(missing(xname) || is.null(xname)) xname <- short.deparse(substitute(x))
  a <- do.call.matched(as.im,
                       list(X=x, ...),
                       c("X", "W",
		         "dimyx", "eps", "xy",
   		         "na.replace", "strict"),
		       sieve=TRUE)
  Z <- a$result
  do.call(hist.im, append(list(x=Z, xname=xname), a$otherargs))
}
