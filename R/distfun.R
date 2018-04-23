#
#   distfun.R
#
#   distance function (returns a function of x,y)
#
#   $Revision: 1.27 $   $Date: 2018/04/23 05:12:30 $
#

distfun <- function(X, ...) {
  UseMethod("distfun")
}

distfun.ppp <- function(X, ..., k=1, undef=Inf) {
  # this line forces X to be bound
  stopifnot(is.ppp(X))
  stopifnot(length(k) == 1)
  force(undef)
  g <- function(x,y=NULL) {
    Y <- xy.coords(x, y)[c("x", "y")]
    if(npoints(X) < k) rep(undef, length(Y$x)) else
    nncross(Y, X, what="dist", k=k) 
  }
  attr(g, "Xclass") <- "ppp"
  g <- funxy(g, as.rectangle(as.owin(X)))
  attr(g, "k") <- k
  attr(g, "extrargs") <- list(k=k, undef=undef)
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
  attr(g, "extrargs") <- list()
  return(g)
}

distfun.owin <- function(X, ..., invert=FALSE) {
  # this line forces X to be bound
  stopifnot(is.owin(X))
  force(invert)
  #
  P <- edges(X)
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
  attr(g, "extrargs") <- list(invert=invert)
  class(g) <- c("distfun", class(g))
  return(g)
}

as.owin.distfun <- function(W, ..., fatal=TRUE) {
  X <- get("X", envir=environment(W))
  result <- if(is.owin(X)) as.rectangle(X) else as.owin(X, ..., fatal=fatal)
  return(result)
}

domain.distfun <- Window.distfun <- function(X, ...) { as.owin(X) }

as.im.distfun <- function(X, W=NULL, ...,
                           eps=NULL, dimyx=NULL, xy=NULL,
                           na.replace=NULL, approx=TRUE) {
  k <- attr(X, "k")
  if(approx && is.null(W) && (is.null(k) || (k == 1))) {
    # use 'distmap' for speed
    env <- environment(X)
    Xdata  <- get("X",      envir=env)
    args <- list(X=Xdata, eps=eps, dimyx=dimyx, xy=xy)
    if(is.owin(Xdata)) {
      args <- append(args, list(invert = get("invert", envir=env)))
    }
    D <- do.call(distmap, args = args)
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
  splat("Distance function for", typestring)
  X <- get("X", envir=environment(x))
  print(X)
  if(!is.null(k <- attr(x, "k")) && k > 1)
    splat("Distance to", ordinal(k), "nearest", objname, "will be computed")
  return(invisible(NULL))
}

summary.distfun <- function(object, ...) {
  xtype <- attr(object, "Xclass")
  w <- as.owin(object)
  fundef <- attr(object, "f")
  attr(fundef, "Xclass") <- NULL
  X <- get("X", envir=environment(object))
  z <- list(xtype   = xtype,
            k       = attr(object, "k") %orifnull% 1,
            Xsumry  = summary(X),
            values  = summary(as.im(object, ...)),
            wintype = w$type,
            frame   = Frame(w),
            units   = unitname(w))
  class(z) <- "summary.distfun"
  return(z)
}

print.summary.distfun <- function(x, ...) {
  typestring <- switch(x$xtype,
                       ppp="point pattern",
                       psp="line segment pattern",
                       owin="window",
                       "unrecognised object")
  objname <- switch(x$xtype,
                    ppp="point",
                    psp="line segment",
                    "object")
  splat("Distance function for", typestring)
  if(x$k > 1)
    splat("Distance to", ordinal(x$k), "nearest", objname, "will be computed")
  windesc <- switch(x$wintype,
                    rectangle="the rectangle",
                    polygonal="a polygonal window inside the frame",
                    mask="a binary mask in the rectangle")
  unitinfo <- summary(x$units)
  sigdig <- getOption('digits')
  splat("defined in",
        windesc,
        prange(signif(x$frame$xrange, sigdig)),
        "x",
        prange(signif(x$frame$yrange, sigdig)),
        unitinfo$plural,
        unitinfo$explain
        )
  v <- x$values
  splat("\nDistance function values:")
  splat("\trange =", prange(signif(v$range, sigdig)))
#  splat("\tintegral =", signif(v$integral, sigdig))
  splat("\tmean =", signif(v$mean, sigdig))
  invisible(NULL)
}

shift.distfun <- rotate.distfun <- scalardilate.distfun <- 
  affine.distfun <- function(X, ...) {
    f <- X
    extrargs <- attr(f, "extrargs")
    if(is.null(extrargs))
      stop(paste("distfun object has outdated format;",
                 "cannot apply geometrical transformation"),
           call.=FALSE)
    Y <- get("X", envir=environment(f))
    Ynew <- do.call(.Generic, list(Y, ...))
    fnew <- do.call(distfun, append(list(Ynew), extrargs))
    return(fnew)
  }


flipxy.distfun <- reflect.distfun <- function(X) {
  f <- X
  extrargs <- attr(f, "extrargs")
  if(is.null(extrargs))
    stop(paste("distfun object has outdated format;",
               "cannot apply geometrical transformation"),
         call.=FALSE)
  Y <- get("X", envir=environment(f))
  Ynew <- do.call(.Generic, list(Y))
  fnew <- do.call(distfun, append(list(Ynew), extrargs))
  return(fnew)
}

rescale.distfun <- function(X, s, unitname) {
  if(missing(s)) s <- NULL
  if(missing(unitname)) unitname <- NULL
  f <- X
  Y <- get("X", envir=environment(f))
  Ynew <- rescale(Y, s, unitname)
  extrargs <- attr(f, "extrargs")
  if(is.null(extrargs))
    stop(paste("distfun object has outdated format;",
               "cannot rescale it"),
         call.=FALSE)
  fnew <- do.call(distfun, append(list(Ynew), extrargs))
  return(fnew)
}
  
