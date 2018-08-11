##
## persp.im.R
##
##  'persp' method for image objects
##      plus annotation
##  
##  $Revision: 1.20 $ $Date: 2016/09/01 05:49:42 $
##

persp.im <- local({

  persp.im <- function(x, ...,
                       colmap=NULL, colin=x, apron=FALSE,
                       visible=FALSE) {
    xname <- deparse(substitute(x))
    xinfo <- summary(x)
    if(xinfo$type == "factor")
      stop("Perspective plot is inappropriate for factor-valued image")
    ## check whether 'col' was specified when 'colmap' was intended
    Col <- list(...)$col
    if(is.null(colmap) && !is.null(Col) && !is.matrix(Col) && length(Col) != 1)
      warning("Argument col is not a matrix. Did you mean colmap?")
    if(!missing(colin)) {
      ## separate image to determine colours
      verifyclass(colin, "im")
      if(!compatible(colin, x)) {
        ## resample 'colin' onto grid of 'x'
        colin <- as.im(colin, W=x)
      }
      if(is.null(colmap))
        colmap <- spatstat.options("image.colfun")(128)
    }
    pop <- spatstat.options("par.persp")
    ##
    if(is.function(colmap) && !inherits(colmap, "colourmap")) {
      ## coerce to a 'colourmap' if possible
      clim <- range(colin, finite=TRUE)
      if(names(formals(colmap))[1] == "n") {
        colval <- colmap(128)
        colmap <- colourmap(colval, range=clim)
      } else {
        ## colour map determined by a rule (e.g. 'beachcolours')
        colmap <- invokeColourmapRule(colmap, colin,
                                      zlim=clim, colargs=list(...))
        if(is.null(colmap))
          stop("Unrecognised syntax for colour function")
      }
    }
    ## colour map?
    if(is.null(colmap)) {
      colinfo <- list(col=NULL)
    } else if(inherits(colmap, "colourmap")) {
      ## colour map object
      ## apply colour function to image data
      colval <- eval.im(colmap(colin))
      colval <- t(as.matrix(colval))
      ## strip one row and column for input to persp.default
      colval <- colval[-1, -1]
      ## replace NA by arbitrary value
      isna <- is.na(colval)
      if(any(isna)) {
        stuff <- attr(colmap, "stuff")
        colvalues <- stuff$outputs
        colval[isna] <- colvalues[1]
      }
      ## pass colour matrix (and suppress lines)
      colinfo <- list(col=colval, border=NA)
    } else {
      ## interpret 'colmap' as colour map
      if(is.list(colmap) && all(c("breaks", "col") %in% names(colmap))) {
        breaks <- colmap$breaks
        colvalues <- colmap$col
      } else if(is.vector(colmap)) {
        colvalues <- colmap
        breaks <- quantile(colin,
                           seq(from=0,to=1,length.out=length(colvalues)+1))
        if(!all(ok <- !duplicated(breaks))) {
          breaks <- breaks[ok]
          colvalues <- colvalues[ok[-1]]
        }
      } else warning("Unrecognised format for colour map")
      ## apply colour map to image values
      colid <- cut.im(colin, breaks=breaks, include.lowest=TRUE)
      colval <- eval.im(colvalues[unclass(colid)])
      colval <- t(as.matrix(colval))
#      nr <- nrow(colval)
#      nc <- ncol(colval)
      ## strip one row and column for input to persp.default
      colval <- colval[-1, -1]
      colval[is.na(colval)] <- colvalues[1]
      ## pass colour matrix (and suppress lines)
      colinfo <- list(col=colval, border=NA)
    }

    if(apron) {
      ## add an 'apron'
      zlim <- list(...)$zlim
      bottom <- if(!is.null(zlim)) zlim[1] else min(x)
      x <- na.handle.im(x, na.replace=bottom)
      x <- padimage(x, bottom)
      xinfo <- summary(x)
      if(is.matrix(colval <- colinfo$col)) {
        colval <- matrix(col2hex(colval), nrow(colval), ncol(colval))
        grijs <- col2hex("lightgrey")
        colval <- cbind(grijs, rbind(grijs, colval, grijs), grijs)
        colinfo$col <- colval
      }
    }

    if(spatstat.options("monochrome"))
      colinfo$col <- to.grey(colinfo$col)
  
    ## get reasonable z scale while fixing x:y aspect ratio
    if(xinfo$type %in% c("integer", "real")) {
      zrange <- xinfo$range
      if(diff(zrange) > 0) {
        xbox <- as.rectangle(x)
        zscale <- 0.5 * mean(diff(xbox$xrange), diff(xbox$yrange))/diff(zrange)
        zlim <- zrange
      } else {
        zscale <- NULL
        mx <- xinfo$mean
        zlim <- mx + c(-1,1) * if(mx == 0) 0.1 else min(abs(mx), 1)
      }
    } else 
      zscale <- zlim <- NULL

    dotargs <- list(...)
    if(spatstat.options("monochrome"))
      dotargs <- col.args.to.grey(dotargs)
    
    yargh <- resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                              dotargs,
                              pop,
                              colinfo,
                              list(xlab="x", ylab="y", zlab=xname),
                              list(scale=FALSE, expand=zscale,
                                   zlim=zlim),
                              list(main=xname),
                              .StripNull=TRUE)

    jawab <- do.call.matched(persp, yargh, 
                             funargs=graphicsPars("persp"))

    attr(jawab, "expand") <- yargh$expand
    
    if(visible)
      attr(jawab, "visible") <- perspvis(x, M=jawab)
    
    return(invisible(jawab))
  }

  diffit <- function(x) {
    y <- diff(x)
    return(c(y[1], y))
  }
  
  perspvis <- function(X, ..., M=NULL) {
    stopifnot(is.im(X))
    ## determine perspective matrix
    if(is.null(M))
      M <- persp(X, ...)
    ## project the coordinates
    ## onto (x,y) plane of plot and z axis pointing out of it
    xy <- rasterxy.im(X, drop=TRUE)
    z <- X[drop=TRUE]
    xyz <- cbind(xy, z)
    v <- cbind(xyz, 1) %*% M
    pxyz <- v[,1:3]/v[,4]
    px <- pxyz[,1]
    py <- pxyz[,2]
    pz <- pxyz[,3]
    ## determine greatest possible difference in 'depth' in one pixel step
    PZ <- as.matrix(X)
    ok <- !is.na(PZ)
    PZ[ok] <- pz
    maxslip <- max(0, abs(apply(PZ, 1, diff)),
                      abs(apply(PZ, 2, diff)), na.rm=TRUE)
    ## determine which pixels are in front
    d <- ceiling(dim(X)/2)
    jx <- cut(px, breaks=d[2])
    iy <- cut(py, breaks=d[1])
    zmax <- tapply(pz, list(iy,jx), max)
    isvis <- infront <- (pz > zmax[cbind(iy,jx)] - 2 * maxslip)
    ##
    if(TRUE) {
      ## Additionally check whether unit normal to surface is pointing to viewer
      Xmat <- as.matrix(X)
      dzdx <- cbind(0, t(apply(Xmat, 1, diff)))/X$xstep
      dzdy <- rbind(0, apply(Xmat, 2, diff))/X$ystep
      dzdx <- as.vector(dzdx[ok])
      dzdy <- as.vector(dzdy[ok])
      ## unscaled normal is (-dzdx, -dzdy, 1)
      if(FALSE) {
        ## THIS DOESN'T WORK - not sure why.
        ## rescale so that length is half diameter of pixel
        fac <- sqrt(X$xstep^2 + X$ystep^2)/(2 * sqrt(dzdx^2+dzdy^2+1))
        ## add to spatial coordinates
        xyzplus <- xyz + fac * cbind(-dzdx, -dzdy, 1)
        ## transform
        vplus <- cbind(xyzplus, 1) %*% M
        pplus <- vplus[,1:3]/vplus[,4]
        ## determine whether normal is pointing toward viewer
        deltaz <- pplus[,3] - pz
        isvis <- infront & (deltaz > 0)
      } else {
        theta <- atan2(M[2,1],M[1,1]) + pi/2
        phi <-  - atan2(M[3,3], M[3,2])
        ## check agreement
        ## cat(paste("Guess: theta=", theta * 180/pi, "\n"))
        ## cat(paste("Guess: phi=", phi * 180/pi, "\n"))
        ## view vector
        viewer <- cos(phi) * c(cos(theta), sin(theta), 0)
                       + c(0, 0, sin(phi))
        ## inner product
        dotprod <- -dzdx * viewer[1] - dzdy * viewer[2] + viewer[3]
        isvis <- infront & (dotprod < 0)
      }
    }
    ## put into image
    Y <- eval.im(X > 0)
    Y[] <- isvis
    ## replace 'NA' by 'FALSE'
    if(anyNA(Y))
      Y <- as.im(Y, na.replace=FALSE)
    return(Y)
  }

  persp.im
})


perspPoints <- function(x, y=NULL, ..., Z, M) {
  xy <- xy.coords(x, y)
  stopifnot(is.im(Z))
  X <- as.ppp(xy, W=Frame(Z))
  if(!(is.matrix(M) && all(dim(M) == 4)))
    stop("M should be a 4 x 4 matrix, returned from persp()")
  V <- attr(M, "visible")
  if(is.null(V)) {
    warning(paste("M does not contain visibility information;",
               "it should be recomputed by persp() with visible=TRUE"))
  } else {
    ## restrict to visible points
    VX <- V[X, drop=FALSE]
    VX[is.na(VX)] <- FALSE
    X <- X[VX]
  }
  #' determine heights
  ZX <- Z[X, drop=FALSE] # may contain NA
  #' transform and plot
  points(trans3d(X$x, X$y, ZX, M), ...)
}

perspSegments <- local({
  perspSegments <- function(x0, y0=NULL, x1=NULL, y1=NULL, ..., Z, M) {
    stopifnot(is.im(Z))
    if(!(is.matrix(M) && all(dim(M) == 4)))
      stop("M should be a 4 x 4 matrix, returned from persp()")
    V <- attr(M, "visible")
    if(is.null(V))
      warning(paste("M does not contain visibility information;",
                 "it should be recomputed by persp() with visible=TRUE"))
    
    if(is.psp(X <- x0) && is.null(y0) && is.null(x1) && is.null(y1)) {
      eX <- X$ends
#      nX <- nrow(eX)
    } else {
#      nX <- length(x0)
      check.nvector(x0, naok=TRUE)
      check.nvector(y0, naok=TRUE)
      check.nvector(x1, naok=TRUE)
      check.nvector(y1, naok=TRUE)
      eX <- cbind(x0, y0, x1, y1)
    }
    if(is.null(V)) {
      Y <- eX
    } else {
      ## chop segments to length of single pixel
      eps <- with(Z, min(xstep,ystep))
      Y <- do.call(rbind, lapply(as.data.frame(t(eX)), chopsegment, eps=eps))
      ## determine which segments are visible
      yleft  <- list(x=Y[,1], y=Y[,2])
      yright <- list(x=Y[,3], y=Y[,4])
      ok <- V[yleft, drop=FALSE] & V[yright, drop=FALSE]
      ok[is.na(ok)] <- FALSE
      Y <- Y[ok, ,drop=FALSE]
    }
    if(nrow(Y) == 0) return(invisible(NULL))
    ## map to projected plane
    x0y0 <- trans3d(Y[,1], Y[,2], Z[list(x=Y[,1],y=Y[,2]), drop=FALSE], M)
    x1y1 <- trans3d(Y[,3], Y[,4], Z[list(x=Y[,3],y=Y[,4]), drop=FALSE], M)
    segments(x0y0$x, x0y0$y, x1y1$x, x1y1$y, ...)
  }

  chopsegment <- function(x, eps) {
    len2 <- (x[3] - x[1])^2 + (x[4] - x[2])^2
    if(len2 <= eps^2) return(x)
    n <- ceiling(sqrt(len2)/eps)
    b <- (1:n)/n
    a <- (0:(n-1))/n
    return(cbind(x[1] + a * (x[3]-x[1]),
                 x[2] + a * (x[4]-x[2]),
                 x[1] + b * (x[3]-x[1]),
                 x[2] + b * (x[4]-x[2])))
  }
      
  perspSegments
})

perspLines <- function(x, y=NULL, ..., Z, M) {
  xy <- xy.coords(x, y)
  n <- length(xy$x)
  perspSegments(x[-n], y[-n], x[-1], y[-1], Z=Z, M=M, ...)
}

perspContour <- function(Z, M, ...,
                         nlevels=10, levels=pretty(range(Z), nlevels)) {
  cl <- contourLines(x=Z$xcol,
                     y=Z$yrow,
                     z=t(Z$v),
                     nlevels=nlevels, levels=levels)
  for(i in seq_along(cl)) {
    cli <- cl[[i]]
    perspLines(cli$x, cli$y, ..., Z=Z, M=M)
  }
  invisible(NULL)
}

