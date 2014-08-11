#
#           pixellate.R
#
#           $Revision: 1.5 $    $Date: 2011/05/18 08:19:08 $
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

pixellate.ppp <- function(x, W=NULL, ..., weights=NULL, padzero=FALSE)
{
    verifyclass(x, "ppp")

    if(!is.null(W))
      W <- as.mask(W)
    else {
      # determine W using as.mask
      dotargs <- list(...)
      namesargs <- names(dotargs)
      matched <- namesargs %in% names(formals(as.mask))
      W <- do.call("as.mask", append(list(x$window), dotargs[matched]))
    } 

    if(x$n == 0) {
      zeroimage <- as.im(as.double(0), W)
      if(padzero) # map NA to 0
        zeroimage <- na.handle.im(zeroimage, 0)
      return(zeroimage)
    }
    
    pixels <- nearest.raster.point(x$x, x$y, W)
    nr <- W$dim[1]
    nc <- W$dim[2]
    if(is.null(weights)) {
    ta <- table(row = factor(pixels$row, levels = 1:nr), col = factor(pixels$col,
        levels = 1:nc))
    } else {
        ta <- tapply(weights, list(row = factor(pixels$row, levels = 1:nr),
                    col = factor(pixels$col, levels=1:nc)), sum)
        ta[is.na(ta)] <- 0
    }
    if(nr == 1 || nc == 1) {
      # 1 x 1 image: need to specify xrange, yrange explicitly
      out <- im(ta, xrange = W$xrange, yrange = W$yrange, unitname=unitname(W))
    } else {
      # normal case: use exact xcol, yrow values
      out <- im(ta, xcol = W$xcol, yrow = W$yrow, unitname=unitname(W))
    }
    # clip to window of data
    if(!padzero)
      out <- out[W, drop=FALSE]
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
  W <- as.mask(W, ...)
  #
  x0 <- W$xrange[1]
  y0 <- W$yrange[1]
  dx <- W$xstep
  dy <- W$ystep
  nx <- W$dim[2]
  ny <- W$dim[1]
  # set up output image (real-valued) and initialise to zero
  Z <- as.im(W, value=pi, na.replace=pi)
  Z <- eval.im(Z * 0)
  # process each component polygon  
  B <- P$bdry
  DUP <- spatstat.options("dupC")
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
             status=as.integer(integer(1)),
             DUP=DUP,
             PACKAGE="spatstat")
    if(zz$status != 0)
      stop("Internal error")
    # increment output image
    Z$v <- Z$v + matrix(zz$out, ny, nx, byrow=TRUE)
  }
  # revert to original scale
  pixelarea <- dx * dy
  Z <- eval.im(Z * pixelarea)
  return(Z)
}

    
  
