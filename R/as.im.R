#
#    as.im.R
#
#    conversion to class "im"
#
#    $Revision: 1.48 $   $Date: 2017/09/04 05:03:30 $
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
  nopar <- is.null(eps) && is.null(dimyx) && is.null(xy)
  if(is.null(W)) {
    if(nopar) {
      X <- repair.image.xycoords(X)
      X <- na.handle.im(X, na.replace)
      return(X)
    }
    # pixel raster determined by dimyx etc
    W <- as.mask(as.rectangle(X), eps=eps, dimyx=dimyx, xy=xy)
    # invoke as.im.owin
    Y <- as.im(W)
  } else if(is.mask(W) || is.im(W) || !nopar) {
    #' raster information is present in { W, eps, dimyx, xy }
    Y <- as.im(W, eps=eps, dimyx=dimyx, xy=xy)
  } else {
    #' use existing raster information in X
    return(X[W, drop=FALSE, tight=TRUE])
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
    v <- matrix(value, d[1L], d[2L])
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
    Z <- im(matrix(0, d[1L], d[2L]), W$xcol, W$yrow, unitname=unitname(X))    
    # adjust values to indicator of X
    Z[X] <- 1
    if(missing(value) && is.null(na.replace)) {
      # done
      out <- Z
    } else {
      # map {0, 1} to {na.replace, value}
      v <- matrix(ifelseAB(Z$v == 0, na.replace, value), d[1L], d[2L])
      out <- im(v, W$xcol, W$yrow, unitname=unitname(X))
    }
    return(out)
  }
  if(X$type == "mask") {
    # raster dimensions determined by X
    # convert X to image
    d <- X$dim
    v <- matrix(value, d[1L], d[2L])
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
  v <- matrix(value, d[1L], d[2L])
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
                           na.replace=NULL, strict=FALSE) {
  f <- X
  if(is.null(W))
    stop("A window W is required")
  W <- as.owin(W)
  W <- as.mask(W, eps=eps, dimyx=dimyx, xy=xy)
  m <- W$m
  funnywindow <- !all(m)

  xx <- as.vector(rasterx.mask(W))
  yy <- as.vector(rastery.mask(W))

  argh <- list(...)
  if(strict) argh <- argh[names(argh) %in% names(formals(f))]

  # evaluate function value at each pixel 
  if(!funnywindow) 
    values <- do.call(f, append(list(xx, yy), argh))
  else {
    # evaluate only inside window
    inside <- as.vector(m)
    val <- do.call(f, append(list(xx[inside], yy[inside]), argh))
    # create space for full matrix
    msize <- length(m)
    values <-
      if(!is.factor(val))
        vector(mode=typeof(val), length=msize)
      else {
        lev <- levels(val)
        factor(rep.int(lev[1L], msize), levels=lev)
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
    return(as.im(X, W, ..., eps=eps, dimyx=dimyx, xy=xy, na.replace=na.replace))
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
      x <- (x[-1L] + x[-lx])/2
    else if(lx != nr)
      stop("length of x coordinate vector does not match number of rows of z")
    if(ly == nc + 1)
      y <- (y[-1L] + y[-ly])/2
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

as.im.data.frame <- function(X, ..., step, fatal=TRUE, drop=TRUE) {
  if(missing(step)) {
    xstep <- ystep <- NULL
  } else {
    step <- ensure2vector(step)
    xstep <- step[1L]
    ystep <- step[2L]
  }
  if(ncol(X) < 3) {
    whinge <- "Argument 'X' must have at least 3 columns of data"
    if(fatal) stop(whinge)
    warning(whinge)
    return(NULL)
  }
  ## extract (x,y) coordinates
  mch <- matchNameOrPosition(c("x", "y", "z"), names(X))
  x <- X[, mch[1L]]
  y <- X[, mch[2L]]
  z <- X[, -mch[1:2], drop=FALSE]
  ## unique x,y coordinates
  xx <- sort(unique(x))
  yy <- sort(unique(y))
  jj <- match(x, xx)
  ii <- match(y, yy)
  iijj <- cbind(ii, jj)
  ## make matrix (for incomplete x, y sequence)
  ok <- checkbigmatrix(length(xx), length(yy), fatal=fatal)
  if(!ok) return(NULL)
  mm <- matrix(NA, length(yy), length(xx))
  ## ensure xx and yy are complete equally-spaced sequences
  fx <- fillseq(xx, step=xstep)
  fy <- fillseq(yy, step=ystep)
  xcol <- fx[[1L]]
  yrow <- fy[[1L]]
  ## trap very large matrices
  ok <- checkbigmatrix(length(xcol), length(yrow), fatal=fatal)
  if(!ok) return(NULL)
  ## mapping from xx to xcol, yy to yrow
  jjj <- fx[[2L]]
  iii <- fy[[2L]]
  ## make matrix for full sequence
  m <- matrix(NA, length(yrow), length(xcol))
  ## run through columns of pixel values
  nz <- ncol(z)
  result <- vector(mode="list", length=nz)
  names(result) <- colnames(z)
  for(k in seq_len(nz)) {
    zk <- z[,k]
    mm[iijj] <- zk
    m[iii,jjj] <- mm
    lev <- levels(zk)
    mo <- if(is.null(lev)) m else factor(as.vector(m), levels=lev)
    result[[k]] <- im(mat=mo, xcol=xcol, yrow=yrow) 
  }
  if(nz == 1 && drop) result <- result[[1L]]
  return(result)
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
  v <- x$v
  if(is.null(dim(v))) 
    dim(v) <- c(length(x$yrow), length(x$xcol))
  im(v, xrange=x$xrange, yrange=x$yrange, unitname=unitname(x))
}
