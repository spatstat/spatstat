#
#           pixellate.R
#
#           $Revision: 1.25 $    $Date: 2017/11/15 07:23:16 $
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

pixellate.ppp <- function(x, W=NULL, ..., weights=NULL, padzero=FALSE,
                          fractional=FALSE, preserve=FALSE,
                          DivideByPixelArea=FALSE) {
  verifyclass(x, "ppp")

  if(is.null(W))
    W <- Window(x)
  isrect <- is.rectangle(W)
  preserve <- preserve && !isrect
  iscount <- is.null(weights) && !fractional && !preserve
  
  W <- do.call.matched(as.mask,
                       resolve.defaults(list(...),
                                        list(w=W)))

  nx <- npoints(x)
  
  insideW <- W$m
  dimW   <- W$dim
  nr <- dimW[1L]
  nc <- dimW[2L]
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
  if(nx == 0) {
    zerovalue <- if(iscount) 0L else as.double(0)
    zeroimage <- as.im(zerovalue, W)
    if(padzero) # map NA to 0
      zeroimage <- na.handle.im(zeroimage, zerovalue)
    result <- zeroimage
    if(k > 1) {
      result <- as.solist(rep(list(zeroimage), k))
      names(result) <- colnames(weights)
    }
    return(result)
  }

  # map points to pixels 
  xx <- x$x
  yy <- x$y
  if(!fractional) {
    #' map (x,y) to nearest raster point
    pixels <- if(preserve) nearest.valid.pixel(xx, yy, W) else 
              nearest.raster.point(xx, yy, W)
    rowfac <- factor(pixels$row, levels=1:nr)
    colfac <- factor(pixels$col, levels=1:nc)
  } else {
    #' attribute fractional weights to the 4 pixel centres surrounding (x,y)
    #' find surrounding pixel centres
    jj <- findInterval(xx, xcolW, rightmost.closed=TRUE)
    ii <- findInterval(yy, yrowW, rightmost.closed=TRUE)
    jleft <- pmax(jj, 1)
    jright <- pmin(jj + 1, nr) 
    ibot <- pmax(ii, 1)
    itop <- pmin(ii+1, nc)
    #' compute fractional weights
    wleft <- pmin(1, abs(xcolW[jright] - xx)/W$xstep)
    wright <- 1 - wleft
    wbot <- pmin(1, abs(yrowW[itop] - yy)/W$ystep)
    wtop <- 1 - wbot
    #' pack together
    ww <- c(wleft * wbot, wleft * wtop, wright * wbot, wright * wtop)
    rowfac <- factor(c(ibot, itop, ibot, itop), levels=1:nr)
    colfac <- factor(c(jleft, jleft, jright, jright), levels=1:nc)
    if(preserve) {
      #' normalise fractions for each data point to sum to 1 inside window
      ok <- insideW[cbind(as.integer(rowfac), as.integer(colfac))]
      wwok <- ww * ok
      denom <- .colSums(wwok, 4, nx, na.rm=TRUE)
      recip <- ifelse(denom == 0, 1, 1/denom)
      ww <- wwok * rep(recip, each=4)
    }
    #' data weights must be replicated
    if(is.null(weights)) {
      weights <- ww
    } else if(k == 1) {
      weights <- ww * rep(weights, 4)
    } else {
      weights <- ww * apply(weights, 2, rep, times=4)
    }
  }
  
  #' sum weights
  if(is.null(weights)) {
    ta <- table(row = rowfac, col = colfac)
  } else if(k == 1) {
    ta <- tapplysum(weights, list(row = rowfac, col=colfac))
  } else {
    ta <- list()
    for(j in 1:k) {
      ta[[j]] <- tapplysum(weights[,j], list(row = rowfac, col=colfac))
    }
  }

  #' divide by pixel area?
  if(DivideByPixelArea) {
    pixelarea <- W$xstep * W$ystep
    if(k == 1) {
      ta <- ta/pixelarea
    } else {
      ta <- lapply(ta, "/", e2=pixelarea)
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
    template <- im(ta[[1L]],
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

pixellate.owin <- function(x, W=NULL, ..., DivideByPixelArea=FALSE) {
  stopifnot(is.owin(x))
  P <- as.polygonal(x)
  R <- as.rectangle(x)
  if(is.null(W)) 
    W <- R
  else if(!is.subset.owin(R, as.rectangle(W)))
    stop("W does not cover the domain of x")
  
  W <- do.call.matched(as.mask,
                       resolve.defaults(list(...),
                                        list(w=W)))
  ## compute
  Zmat <- polytileareaEngine(P, W$xrange, W$yrange, nx=W$dim[2L], ny=W$dim[1L],
                             DivideByPixelArea)
  ## convert to image
  Z <- im(Zmat, xcol=W$xcol, yrow=W$yrow, xrange=W$xrange, yrange=W$yrange,
          unitname=unitname(W))
  return(Z)
}

polytileareaEngine <- function(P, xrange, yrange, nx, ny,
                               DivideByPixelArea=FALSE) {
  x0 <- xrange[1L]
  y0 <- yrange[1L]
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
    xx <- c(xx, xx[1L])
    yy <- c(yy, yy[1L])
    nn <- nn+1
    # call C routine
    zz <- .C("poly2imA",
             ncol=as.integer(nx),
             nrow=as.integer(ny),
             xpoly=as.double(xx),
             ypoly=as.double(yy),
             npoly=as.integer(nn),
             out=as.double(numeric(nx * ny)),
             status=as.integer(integer(1L)),
             PACKAGE = "spatstat")
    if(zz$status != 0)
      stop("Internal error")
    # increment output 
    Z[] <- Z[] + zz$out
  }
  if(!DivideByPixelArea) {
    #' revert to original scale
    pixelarea <- dx * dy
    Z <- Z * pixelarea
  }
  return(Z)
}


  
