#
#    as.im.R
#
#    conversion to class "im"
#
#    $Revision: 1.42 $   $Date: 2015/04/27 08:51:04 $
#
#    as.im()
#

as.im <- function(X, ...) {
  UseMethod("as.im")
}

as.im.im <- function(X, W=NULL, ...,
                     eps=NULL, dimyx=NULL, xy=NULL,
                     na.replace=NULL) {
  X <- repair.old.factor.image(X)
  if(is.null(W)) {
    if(is.null(eps) && is.null(dimyx) && is.null(xy)) {
      X <- repair.image.xycoords(X)
      X <- na.handle.im(X, na.replace)
      return(X)
    }
    # pixel raster determined by dimyx etc
    W <- as.mask(as.rectangle(X), eps=eps, dimyx=dimyx, xy=xy)
    # invoke as.im.owin
    Y <- as.im(W)
  } else {
    # apply dimyx (etc) if present,
    # otherwise use W to determine pixel raster
    Y <- as.im(W, eps=eps, dimyx=dimyx, xy=xy)
  }
  # resample X onto raster of Y
  Y <- rastersample(X, Y)

  return(na.handle.im(Y, na.replace))
}

as.im.owin <- function(X, W=NULL, ...,
                       eps=NULL, dimyx=NULL, xy=NULL,
                       na.replace=NULL, value=1) {
  if(!(is.null(eps) && is.null(dimyx) && is.null(xy))) {
    # raster dimensions determined by dimyx etc
    # convert X to a mask 
    M <- as.mask(X, eps=eps, dimyx=dimyx, xy=xy)
    # convert mask to image
    d <- M$dim
    v <- matrix(value, d[1], d[2])
    m <- M$m
    v[!m] <- if(is.null(na.replace)) NA else na.replace
    out <- im(v, M$xcol, M$yrow,
              xrange=M$xrange, yrange=M$yrange,
              unitname=unitname(X))
    return(out)
  }
  if(!is.null(W) && is.owin(W) && W$type == "mask") {
    # raster dimensions determined by W
    # convert W to zero image
    d <- W$dim
    Z <- im(matrix(0, d[1], d[2]), W$xcol, W$yrow, unitname=unitname(X))    
    # adjust values to indicator of X
    Z[X] <- 1
    if(missing(value) && is.null(na.replace)) {
      # done
      out <- Z
    } else {
      # map {0, 1} to {na.replace, value}
      v <- matrix(ifelseAB(Z$v == 0, na.replace, value), d[1], d[2])
      out <- im(v, W$xcol, W$yrow, unitname=unitname(X))
    }
    return(out)
  }
  if(X$type == "mask") {
    # raster dimensions determined by X
    # convert X to image
    d <- X$dim
    v <- matrix(value, d[1], d[2])
    m <- X$m
    v[!m] <- if(is.null(na.replace)) NA else na.replace
    out <- im(v, xcol=X$xcol, yrow=X$yrow,
              xrange=X$xrange, yrange=X$yrange, unitname=unitname(X))
    return(out)
  }
  # X is not a mask.
  # W is either missing, or is not a mask.
  # Convert X to a image using default settings
  M <- as.mask(X)
  # convert mask to image
  d <- M$dim
  v <- matrix(value, d[1], d[2])
  m <- M$m
  v[!m] <- if(is.null(na.replace)) NA else na.replace
  out <- im(v, M$xcol, M$yrow, unitname=unitname(X))
  return(out)
}

as.im.funxy <- function(X, W=Window(X), ...) {
  as.im.function(X, W=W, ...)
}

as.im.function <- function(X, W=NULL, ...,
                           eps=NULL, dimyx=NULL, xy=NULL,
                           na.replace=NULL) {
  f <- X
  if(is.null(W))
    stop("A window W is required")
  W <- as.owin(W)
  W <- as.mask(W, eps=eps, dimyx=dimyx, xy=xy)
  m <- W$m
  funnywindow <- !all(m)

  xx <- as.vector(rasterx.mask(W))
  yy <- as.vector(rastery.mask(W))

  # evaluate function value at each pixel 
  if(!funnywindow) 
    values <- f(xx, yy, ...)
  else {
    # evaluate only inside window
    inside <- as.vector(m)
    val <- f(xx[inside], yy[inside], ...)
    # create space for full matrix
    msize <- length(m)
    values <-
      if(!is.factor(val))
        vector(mode=typeof(val), length=msize)
      else {
        lev <- levels(val)
        factor(rep.int(lev[1], msize), levels=lev)
      }
    # copy values, assigning NA outside window
    values[inside] <- val
    values[!inside] <- NA
  }

  nc <- length(W$xcol)
  nr <- length(W$yrow)
  if(nr == 1 || nc == 1) {
    # exception: can't determine pixel width/height from centres
    out <- im(matrix(values, nr, nc),
              xrange=W$xrange, yrange=W$yrange, unitname=unitname(W))
  } else {
    out <- im(values, W$xcol, W$yrow, unitname=unitname(W))
  }
  return(na.handle.im(out, na.replace))
}

as.im.matrix <- function(X, W=NULL, ...) {
  nr <- nrow(X)
  nc <- ncol(X)
  if(is.null(W))
    return(im(X, ...))
  W <- as.owin(W)
  if(W$type == "mask") {
    xcol <- W$xcol
    yrow <- W$yrow
    # pixel coordinate information
    if(length(xcol) == nc && length(yrow) == nr)
      return(im(X, xcol, yrow, unitname=unitname(W)))
  }
  # range information
  R <- as.rectangle(W)
  xrange <- R$xrange
  yrange <- R$yrange
  return(im(X, xrange=xrange, yrange=yrange, unitname=unitname(W)))
}

as.im.default <- function(X, W=NULL, ...,
                          eps=NULL, dimyx=NULL, xy=NULL,
                          na.replace=NULL) {

  if((is.vector(X) || is.factor(X)) && length(X) == 1) {
    # numerical value: interpret as constant function
    xvalue <- X
    X <- function(xx, yy, ...) { rep.int(xvalue, length(xx)) }
    return(as.im(X, W, ..., dimyx=dimyx, na.replace=na.replace))
  }
  
  if(is.list(X) && checkfields(X, c("x","y","z"))) {
    stopifnot(is.matrix(X$z))
    z <- X$z
    y <- X$y
    x <- X$x
    # Usual S convention as in contour.default() and image.default()
    # Rows of z correspond to x values.
    nr <- nrow(z)
    nc <- ncol(z)
    lx <- length(x)
    ly <- length(y)
    if(lx == nr + 1)
      x <- (x[-1] + x[-lx])/2
    else if(lx != nr)
      stop("length of x coordinate vector does not match number of rows of z")
    if(ly == nc + 1)
      y <- (y[-1] + y[-ly])/2
    else if(ly != nc)
      stop("length of y coordinate vector does not match number of columns of z")
    # convert to class "im"
    out <- im(t(z), x, y)
    # now apply W and dimyx if present
    if(is.null(W) && !(is.null(eps) && is.null(dimyx) && is.null(xy)))
      out <- as.im(out, eps=eps, dimyx=dimyx, xy=xy)
    else if(!is.null(W))
      out <- as.im(out, W=W, eps=eps, dimyx=dimyx, xy=xy)
    return(na.handle.im(out, na.replace))
  }
  stop("Can't convert X to a pixel image")
}

as.im.ppp <- function(X, ...) {
  pixellate(X, ..., weights=NULL, zeropad=FALSE)
}

# convert to image from some other format, then do something

do.as.im <- function(x, action, ...,
                     W = NULL, eps = NULL, dimyx = NULL, xy = NULL, 
                     na.replace = NULL) {
  Z <- as.im(x, W=W, eps=eps, dimyx=dimyx, xy=xy, na.replace=na.replace)
  Y <- do.call(action, list(Z, ...))
  return(Y)
}

na.handle.im <- function(X, na.replace) {
if(is.null(na.replace))
  return(X)
if(length(na.replace) != 1)
  stop("na.replace should be a single value")
X$v[is.na(X$v)] <- na.replace
return(X)
}

repair.old.factor.image <- function(x) {
  # convert from old to new representation of factor images
  if(x$type != "factor")
    return(x)
  v <- x$v
  isold <- !is.null(lev <- attr(x, "levels"))
  isnew <- is.factor(v) && is.matrix(v)
  if(isnew)
    return(x)
  if(!isold)
    stop("Internal error: unrecognised format for factor-valued image")
  v <- factor(v, levels=lev)
  dim(v) <- x$dim
  x$v <- v
  return(x)
}

repair.image.xycoords <- function(x) {
  im(x$v, xrange=x$xrange, yrange=x$yrange, unitname=unitname(x))
}
