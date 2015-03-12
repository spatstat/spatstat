#
#           pixellate.R
#
#           $Revision: 1.17 $    $Date: 2015/03/12 11:23:56 $
#
#     pixellate            convert an object to a pixel image
#
#     pixellate.ppp        convert a point pattern to a pixel image
#                          (pixel value = number of points in pixel)
#
#     pixellate.owin       convert a window to a pixel image
#                          (pixel value = area of intersection with pixel)
#

pixellate <- function(x, ...) {
  UseMethod("pixellate")
}

pixellate.ppp <- function(x, W=NULL, ..., weights=NULL, padzero=FALSE) {
  verifyclass(x, "ppp")

  if(is.null(W))
    W <- Window(x)

  W <- do.call.matched("as.mask",
                       resolve.defaults(list(...),
                                        list(w=W)))
  
  insideW <- W$m
  dimW   <- W$dim
  xcolW <- W$xcol
  yrowW <- W$yrow
  xrangeW <- W$xrange
  yrangeW <- W$yrange
  unitsW <- unitname(W)
    
  # multiple columns of weights?
  if(is.data.frame(weights) || is.matrix(weights)) {
    k <- ncol(weights)
    stopifnot(nrow(weights) == npoints(x))
    weights <- if(k == 1) as.vector(weights) else as.data.frame(weights)
  } else {
    k <- 1
    if(length(weights) == 0) weights <- NULL else 
      stopifnot(length(weights) == npoints(x) || length(weights) == 1)
    if(length(weights) == 1)
      weights <- rep(weights, npoints(x))
  }

  # handle empty point pattern
  if(x$n == 0) {
    zeroimage <- as.im(as.double(0), W)
    if(padzero) # map NA to 0
      zeroimage <- na.handle.im(zeroimage, 0)
    result <- zeroimage
    if(k > 1) {
      result <- as.solist(rep(list(zeroimage), k))
      names(result) <- colnames(weights)
    }
    return(result)
  }

  # perform calculation
  pixels <- nearest.raster.point(x$x, x$y, W)
  nr <- dimW[1]
  nc <- dimW[2]
  rowfac <- factor(pixels$row, levels=1:nr)
  colfac <- factor(pixels$col, levels=1:nc)
  if(is.null(weights)) {
    ta <- table(row = rowfac, col = colfac)
  } else if(k == 1) {
    ta <- tapply(weights, list(row = rowfac, col=colfac), sum)
    ta[is.na(ta)] <- 0
  } else {
    ta <- list()
    for(j in 1:k) {
      taj <- tapply(weights[,j], list(row = rowfac, col=colfac), sum)
      taj[is.na(taj)] <- 0
      ta[[j]] <- taj
    }
  }

  # pack up as image(s)
  if(k == 1) {
    # single image
    # clip to window of data
    if(!padzero)
      ta[!insideW] <- NA
    out <- im(ta,
              xcol = xcolW, yrow = yrowW,
              xrange = xrangeW, yrange = yrangeW,
              unitname=unitsW)
  } else {
    # case k > 1
    # create template image to reduce overhead
    template <- im(ta[[1]],
                   xcol = xcolW, yrow = yrowW,
                   xrange = xrangeW, yrange = yrangeW,
                   unitname=unitsW)
    out <- list()
    for(j in 1:k) {
      taj <- ta[[j]]
      # clip to window of data
      if(!padzero) 
        taj[!insideW] <- NA
      # copy template and reassign pixel values
      outj <- template
      outj$v <- taj
      # store
      out[[j]] <- outj
    }
    out <- as.solist(out)
    names(out) <- names(weights)
  }
  return(out)
}

pixellate.owin <- function(x, W=NULL, ...) {
  stopifnot(is.owin(x))
  P <- as.polygonal(x)
  R <- as.rectangle(x)
  if(is.null(W)) 
    W <- R
  else if(!is.subset.owin(R, as.rectangle(W)))
    stop("W does not cover the domain of x")
  
  W <- do.call.matched("as.mask",
                       resolve.defaults(list(...),
                                        list(w=W)))
  ## compute
  Zmat <- polytileareaEngine(P, W$xrange, W$yrange, nx=W$dim[2], ny=W$dim[1])
  ## convert to image
  Z <- im(Zmat, xcol=W$xcol, yrow=W$yrow, xrange=W$xrange, yrange=W$yrange,
          unitname=unitname(W))
  return(Z)
}

polytileareaEngine <- function(P, xrange, yrange, nx, ny) {
  x0 <- xrange[1]
  y0 <- yrange[1]
  dx <- diff(xrange)/nx
  dy <- diff(yrange)/ny
  # process each component polygon
  Z <- matrix(0.0, ny, nx)
  B <- P$bdry
  for(i in seq_along(B)) {
    PP <- B[[i]]
    # transform so that pixels become unit squares
    QQ <- affinexypolygon(PP, vec = c(-x0, -y0))
    RR <- affinexypolygon(QQ, mat = diag(1/c(dx, dy)))
    # 
    xx <- RR$x
    yy <- RR$y
    nn <- length(xx)
    # close polygon
    xx <- c(xx, xx[1])
    yy <- c(yy, yy[1])
    nn <- nn+1
    # call C routine
    zz <- .C("poly2imA",
             ncol=as.integer(nx),
             nrow=as.integer(ny),
             xpoly=as.double(xx),
             ypoly=as.double(yy),
             npoly=as.integer(nn),
             out=as.double(numeric(nx * ny)),
             status=as.integer(integer(1)))
    if(zz$status != 0)
      stop("Internal error")
    # increment output 
    Z[] <- Z[] + zz$out
  }
  # revert to original scale
  pixelarea <- dx * dy
  return(Z * pixelarea)
}


  
