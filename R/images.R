#
#       images.R
#
#         $Revision: 1.103 $     $Date: 2013/07/25 09:58:59 $
#
#      The class "im" of raster images
#
#     im()     object creator
#
#     is.im()   tests class membership
#
#     rasterx.im(), rastery.im()    
#                      raster X and Y coordinates
#
#     nearest.pixel()   
#     lookup.im()
#                      facilities for looking up pixel values
#
################################################################
########   basic support for class "im"
################################################################
#
#   creator 

im <- function(mat, xcol=seq_len(ncol(mat)), yrow=seq_len(nrow(mat)), 
               xrange=NULL, yrange=NULL,
               unitname=NULL) {

  typ <- typeof(mat)
  if(typ == "double")
    typ <- "real"

  miss.xcol <- missing(xcol)
  miss.yrow <- missing(yrow)
  
  # determine dimensions
  if(is.matrix(mat)) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    if(length(xcol) != nc)
      stop("Length of xcol does not match ncol(mat)")
    if(length(yrow) != nr)
      stop("Length of yrow does not match nrow(mat)")
  } else {
    if(miss.xcol || miss.yrow)
      stop(paste(sQuote("mat"),
                 "is not a matrix and I can't guess its dimensions"))
    stopifnot(length(mat) == length(xcol) * length(yrow))
    nc <- length(xcol)
    nr <- length(yrow)
  }

  # deal with factor case
  if(is.factor(mat)) {
    typ <- "factor"
  } else if(!is.null(lev <- levels(mat))) {
    typ <- "factor"
    mat <- factor(mat, levels=lev)
  }

  # Ensure 'mat' is a matrix (without destroying factor information)
  if(!is.matrix(mat))
    dim(mat) <- c(nr, nc)

  # set up coordinates
  if((miss.xcol || length(xcol) <= 1) && !is.null(xrange) ) {
    # use 'xrange' 
    xstep <- diff(xrange)/nc
    xcol <- seq(from=xrange[1] + xstep/2, to=xrange[2] - xstep/2, length.out=nc)
  } else if(length(xcol) > 1) {
    # use 'xcol'
    # ensure spacing is constant
    xcol <- seq(from=min(xcol), to=max(xcol), length.out=length(xcol))
    xstep <- diff(xcol)[1]
    xrange <- range(xcol) + c(-1,1) * xstep/2
  } else stop("Cannot determine pixel width")
  
  if((miss.yrow || length(yrow) <= 1) && !is.null(yrange)) {
    # use 'yrange'
    ystep <- diff(yrange)/nr
    yrow <- seq(from=yrange[1] + ystep/2, to=yrange[2] - ystep/2, length.out=nr)
  } else if(length(yrow) > 1) {
    # use 'yrow'
    # ensure spacing is constant
    yrow <- seq(from=min(yrow), to=max(yrow), length.out=length(yrow))
    ystep <- diff(yrow)[1]
    yrange <- range(yrow) + c(-1,1) * ystep/2
  }  else stop("Cannot determine pixel height")

  unitname <- as.units(unitname)

  # get rid of those annoying 8.67e-19 printouts
  swat <- function(x) {ifelseAX(abs(x) < .Machine$double.eps, 0, x)}
  xrange <- swat(xrange)
  yrange <- swat(yrange)
  
  out <- list(v   = mat,
              dim = c(nr, nc),
              xrange   = xrange,
              yrange   = yrange,
              xstep    = xstep,
              ystep    = ystep,
              xcol    = xcol,
              yrow    = yrow,
              type    = typ,
              units   = unitname)
  class(out) <- "im"
  return(out)
}

is.im <- function(x) {
  inherits(x,"im")
}

levels.im <- function(x) {
  levels(x$v)
}

"levels<-.im" <- function(x, value) {
  if(x$type != "factor") 
    stop("image is not factor-valued")
  levels(x$v) <- value
  x
}

################################################################
########   methods for class "im"
################################################################

shift.im <- function(X, vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "im")
  if(!is.null(origin)) {
    stopifnot(is.character(origin))
    if(!missing(vec))
      warning("argument vec ignored; overruled by argument origin")
    origin <- pickoption("origin", origin, c(centroid="centroid",
                                             midpoint="midpoint",
                                             bottomleft="bottomleft"))
    W <- as.owin(X)
    locn <- switch(origin,
                   centroid={ unlist(centroid.owin(W)) },
                   midpoint={ c(mean(W$xrange), mean(W$yrange)) },
                   bottomleft={ c(W$xrange[1], W$yrange[1]) })
    return(shift(X, -locn))
  }
  X$xrange <- X$xrange + vec[1]
  X$yrange <- X$yrange + vec[2]
  X$xcol <- X$xcol + vec[1]
  X$yrow <- X$yrow + vec[2]
  return(X)
}

"[.im" <- 
function(x, i, j, ..., drop=TRUE, raster=NULL, rescue=is.owin(i)) {

  # detect 'blank' arguments like second argument in x[i, ] 
  ngiven <- length(sys.call())
  nmatched <- length(match.call())
  nblank <- ngiven - nmatched
  itype <- if(missing(i)) "missing" else "given"
  jtype <- if(missing(j)) "missing" else "given"
  if(nblank == 1) {
    if(!missing(i)) jtype <- "blank"
    if(!missing(j)) itype <- "blank"
  } else if(nblank == 2) {
    itype <- jtype <- "blank"
  }

  if(missing(rescue) && itype != "given")
    rescue <- FALSE

  if(itype == "missing" && jtype == "missing") {
    # no indices: return entire image 
    out <- if(is.null(raster)) x else as.im(raster)
    xy <- expand.grid(y=out$yrow,x=out$xcol)
    if(!is.null(raster)) {
      # resample image on new pixel raster
      values <- lookup.im(x, xy$x, xy$y, naok=TRUE)
      out <- im(values, out$xcol, out$yrow, unitname=unitname(out))
    }
    if(!drop)
      return(out)
    else {
      v <- out$v
      return(v[!is.na(v)])
    }
  }

  if(itype == "given") {
    # .................................................................
    # Try spatial index
    # .................................................................
    if(verifyclass(i, "owin", fatal=FALSE)) {

      if(jtype == "given")
        warning("Argument j ignored")
      
      # 'i' is a window
      # if drop = FALSE, just set values outside window to NA
      # if drop = TRUE, extract values for all pixels inside window
      #                 as an image (if 'i' is a rectangle)
      #                 or as a vector (otherwise)

      out <- if(is.null(raster)) x else as.im(raster)
      xy <- expand.grid(y=out$yrow,x=out$xcol)
      if(!is.null(raster)) {
        # resample image on new pixel raster
        values <- lookup.im(x, xy$x, xy$y, naok=TRUE)
        out <- im(values, out$xcol, out$yrow, unitname=unitname(out))
      }
      inside <- inside.owin(xy$x, xy$y, i)
      if(!drop) { 
        out$v[!inside] <- NA
        return(out)
      } else if(!rescue || i$type != "rectangle") {
        values <- out$v[inside]
        return(values)
      } else {
        disjoint <- function(r, s) { (r[2] < s[1]) || (r[1] > s[2])  }
        clip <- function(r, s) { c(max(r[1],s[1]), min(r[2],s[2])) }
        inrange <- function(x, r) { (x >= r[1]) & (x <= r[2]) }
        if(disjoint(i$xrange, x$xrange) || disjoint(i$yrange, x$yrange))
          # empty intersection
          return(numeric(0))
        xr <- clip(i$xrange, x$xrange)
        yr <- clip(i$yrange, x$yrange)
        colsub <- inrange(out$xcol, xr)
        rowsub <- inrange(out$yrow, yr)
        ncolsub <- sum(colsub)
        nrowsub <- sum(rowsub)
        if(ncolsub == 0 || nrowsub == 0)
          return(numeric(0))
        marg <- list(mat=out$v[rowsub, colsub, drop=FALSE],
                     unitname=unitname(x))
        xarg <-
          if(ncolsub > 1) list(xcol = out$xcol[colsub]) else list(xrange=xr)
        yarg <-
          if(nrowsub > 1) list(yrow = out$yrow[rowsub]) else list(yrange=yr)
        result <- do.call("im", c(marg, xarg, yarg))
        return(result)
      } 
    }
    if(verifyclass(i, "im", fatal=FALSE)) {
      if(jtype == "given")
        warning("Argument j ignored")
      # logical images OK
      if(i$type == "logical") {
        # convert to window
        w <- as.owin(eval.im(ifelse1NA(i)))
        return(x[w, drop=drop, ..., raster=raster])
      } else stop("Subset argument \'i\' is an image, but not of logical type")
    }

    if(is.ppp(i)) {
      # 'i' is a point pattern 
      if(jtype == "given")
        warning("Argument j ignored")
      # Look up the greyscale values for the points of the pattern
      values <- lookup.im(x, i$x, i$y, naok=TRUE)
      if(drop) 
        values <- values[!is.na(values)]
      if(length(values) == 0) 
        # ensure the zero-length vector is of the right type
        values <- 
          switch(x$type,
                 factor={ factor(, levels=levels(x)) },
                 integer = { integer(0) },
                 logical = { logical(0) },
                 real = { numeric(0) },
                 complex = { complex(0) },
                 character = { character(0) },
                 { values }
                 )
      return(values)
    }
  }
  # ............... not a spatial index .............................

  # Try indexing as a matrix

  # Construct a matrix index call for possible re-use
  M <- as.matrix(x)
  ycall <- switch(itype,
                  given = {
                    switch(jtype,
                           given   = quote(M[i, j, drop=FALSE]),
                           blank   = quote(M[i,  , drop=FALSE]),
                           missing = quote(M[i,    drop=FALSE]))
                  },
                  blank = {
                    switch(jtype,
                           given   = quote(M[ , j, drop=FALSE]),
                           blank   = quote(M[ ,  , drop=FALSE]),
                           missing = quote(M[ ,    drop=FALSE]))
                  },
                  missing = {
                    switch(jtype,
                           given   = quote(M[j=j,  drop=FALSE]),
                           blank   = quote(M[j= ,  drop=FALSE]),
                           missing = quote(M[      drop=FALSE]))
                  })
  # try it
  y <- try(eval(as.call(ycall)), silent=TRUE)
  if(!inherits(y, "try-error")) {
    # valid subset index for a matrix
    if(rescue) {
      # check whether it's a rectangular block, in correct order
      RR <- row(x$v)
      CC <- col(x$v)
      rcall <- ycall
      rcall[[2]] <- quote(RR)
      ccall <- ycall
      ccall[[2]] <- quote(CC)
      rr <- eval(as.call(rcall))
      cc <- eval(as.call(ccall))
      rseq <- sort(unique(as.vector(rr)))
      cseq <- sort(unique(as.vector(cc)))
      if(all(diff(rseq) == 1) && all(diff(cseq) == 1) &&
         (length(rr) == length(rseq) * length(cseq)) &&
         all(rr == RR[rseq, cseq]) && all(cc == CC[rseq,cseq])) {
        # yes - make image
        dim(y) <- c(length(rseq), length(cseq))
        Y <- x
        Y$v <- y
        Y$dim <- dim(y)
        Y$xcol <- x$xcol[cseq]
        Y$yrow <- x$yrow[rseq]
        Y$xrange <- range(Y$xcol) + c(-1,1) * x$xstep/2
        Y$yrange <- range(Y$yrow) + c(-1,1) * x$ystep/2
        return(Y)
      }
    }
    # return pixel values (possibly as matrix)
    return(y)
  }

  # Last chance!
  if(itype == "given" &&
     !is.matrix(i) &&
     !is.null(ip <- as.ppp(i, W=as.owin(x), fatal=FALSE, check=FALSE))) {
    # 'i' is convertible to a point pattern 
    # Look up the greyscale values for the points of the pattern
    values <- lookup.im(x, ip$x, ip$y, naok=TRUE)
    if(drop) 
      values <- values[!is.na(values)]
    if(length(values) == 0) 
      # ensure the zero-length vector is of the right type
      values <- 
        switch(x$type,
               factor={ factor(, levels=levels(x)) },
               integer = { integer(0) },
               logical = { logical(0) },
               real = { numeric(0) },
               complex = { complex(0) },
               character = { character(0) },
               { values }
               )
    return(values)
  }
  
  stop("The subset operation is undefined for this type of index")
}


"[<-.im" <- function(x, i, j, value) {
  
  # detect 'blank' arguments like second argument of x[i, ] 
  ngiven <- length(sys.call())
  nmatched <- length(match.call())
  nblank <- ngiven - nmatched
  itype <- if(missing(i)) "missing" else "given"
  jtype <- if(missing(j)) "missing" else "given"
  if(nblank == 1) {
    if(!missing(i)) jtype <- "blank"
    if(!missing(j)) itype <- "blank"
  } else if(nblank == 2) {
    itype <- jtype <- "blank"
  }

  X <- x
  W <- as.owin(X)
  if(is.im(value)) {
    value <- value$v
  }
  stopifnot(is.vector(value) || is.matrix(value) || is.factor(value))

  if(itype == "missing" && jtype == "missing") {
    # no index provided
    # set all pixels to 'value'
    v <- X$v
    v[!is.na(v)] <- value
    X$v <- v
    return(X)
  }
  if(itype == "given") {
    # ..................... Try a spatial index ....................
    if(verifyclass(i, "owin", fatal=FALSE)) {
      if(jtype == "given") warning("Index j ignored")
      # 'i' is a window
      if(is.empty(i))
        return(X)
      xx <- as.vector(raster.x(W))
      yy <- as.vector(raster.y(W))
      ok <- inside.owin(xx, yy, i)
      X$v[ok] <- value
      return(X)
    }
    if(verifyclass(i, "im", fatal=FALSE) && i$type == "logical") {
      if(jtype == "given") warning("Index j ignored")
      # convert logical vector to window where entries are TRUE
      i <- as.owin(eval.im(ifelse1NA(i)))
      # continue as above
      xx <- as.vector(raster.x(W))
      yy <- as.vector(raster.y(W))
      ok <- inside.owin(xx, yy, i)
      X$v[ok] <- value
      return(X)
    }
    if(is.ppp(i)) {
      # 'i' is a point pattern
      if(jtype == "given") warning("Index j ignored")
      nv <- length(value)
      np <- npoints(i)
      if(nv != np && nv != 1)
        stop("Length of replacement value != number of point locations")
      # test whether all points are inside window FRAME
      ok <- inside.owin(i$x, i$y, as.rectangle(W))
      if(any(!ok)) {
        warning("Some points are outside the outer frame of the image")
        if(nv == np)
          value <- value[ok]
        i <- i[ok]
      }
      if(npoints(i) > 0) {
        # determine row & column positions for each point 
        loc <- nearest.pixel(i$x, i$y, X)
        # set values
        X$v[cbind(loc$row, loc$col)] <- value
      }
      return(X)
    }
  }
  # .................. 'i' is not a spatial index ....................
  
  # Construct a matrix replacement call 
  ycall <- switch(itype,
                  given = {
                    switch(jtype,
                           given   = quote(X$v[i, j] <- value),
                           blank   = quote(X$v[i,  ] <- value),
                           missing = quote(X$v[i]    <- value))
                  },
                  blank = {
                    switch(jtype,
                           given   = quote(X$v[ , j] <- value),
                           blank   = quote(X$v[ ,  ] <- value),
                           missing = quote(X$v[ ] <- value))
                  },
                  missing = {
                    switch(jtype,
                           given   = quote(X$v[j=j] <- value),
                           blank   = quote(X$v[j= ] <- value),
                           missing = quote(X$v[] <- value))
                  })
  # try it
  litmus <- try(eval(as.call(ycall)), silent=TRUE)
  if(!inherits(litmus, "try-error")) 
    return(X)

  #  Last chance!
  if(itype == "given" &&
     !is.matrix(i) &&
     !is.null(ip <- as.ppp(i, W=W, fatal=FALSE, check=TRUE))) {
    # 'i' is convertible to a point pattern
    if(jtype == "given") warning("Index j ignored")
    nv <- length(value)
    np <- npoints(ip)
    if(nv != np && nv != 1)
      stop("Length of replacement value != number of point locations")
    # test whether all points are inside window FRAME
    ok <- inside.owin(ip$x, ip$y, as.rectangle(W))
    if(any(!ok)) {
      warning("Some points are outside the outer frame of the image")
      if(nv == np)
        value <- value[ok]
      ip <- ip[ok]
    }
    if(npoints(ip) > 0) {
      # determine row & column positions for each point 
      loc <- nearest.pixel(ip$x, ip$y, X)
      # set values
      X$v[cbind(loc$row, loc$col)] <- value
    }
    return(X)
  }

  stop("The subset operation is undefined for this type of index")
}

################################################################
########   other tools
################################################################

#
# This function is similar to nearest.raster.point except for
# the third argument 'im' and the different idiom for calculating
# row & column - which could be used in nearest.raster.point()

nearest.pixel <- function(x,y,im) {
  verifyclass(im, "im")
  if(length(x) > 0) {
    nr <- im$dim[1]
    nc <- im$dim[2]
    cc <- round(1 + (x - im$xcol[1])/im$xstep)
    rr <- round(1 + (y - im$yrow[1])/im$ystep)
    cc <- pmax.int(1,pmin.int(cc, nc))
    rr <- pmax.int(1,pmin.int(rr, nr))
  } else cc <- rr <- integer(0)
  return(list(row=rr, col=cc))
}

# Explores the 3 x 3 neighbourhood of nearest.pixel
# and finds the nearest pixel that is not NA

nearest.valid.pixel <- function(x,y,im) {
  rc <- nearest.pixel(x,y,im)
  rr <- rc$row
  cc <- rc$col
  # check whether any pixels are outside image domain
  outside <- is.na(im$v)
  miss <- outside[cbind(rr, cc)]
  if(!any(miss))
    return(rc)
  # for offending pixels, explore 3 x 3 neighbourhood
  nr <- im$dim[1]
  nc <- im$dim[2]
  xcol <- im$xcol
  yrow <- im$yrow
  for(i in which(miss)) {
    rows <- rr[i] + c(-1,0,1)
    cols <- cc[i] + c(-1,0,1)
    rows <- unique(pmax.int(1, pmin.int(rows, nr)))
    cols <- unique(pmax.int(1, pmin.int(cols, nc)))
    rcp <- expand.grid(row=rows, col=cols)
    ok <- !outside[as.matrix(rcp)]
    if(any(ok)) {
      # At least one of the neighbours is valid
      # Find the closest one
      rcp <- rcp[ok,]
      dsq <- with(rcp, (x[i] - xcol[col])^2 + (y[i] - yrow[row])^2)
      j <- which.min(dsq)
      rc$row[i] <- rcp$row[j]
      rc$col[i] <- rcp$col[j]
    }
  }
  return(rc)
}
  

# This function is a generalisation of inside.owin()
# to images other than binary-valued images.

lookup.im <- function(Z, x, y, naok=FALSE, strict=TRUE) {
  verifyclass(Z, "im")

  if(Z$type == "factor")
    Z <- repair.old.factor.image(Z)
  
  if(length(x) != length(y))
    stop("x and y must be numeric vectors of equal length")

  # initialise answer to NA 
  if(Z$type != "factor") {
    niets <- NA
    mode(niets) <- mode(Z$v)
  } else {
    niets <- factor(NA, levels=levels(Z))
  }
  value <- rep.int(niets, length(x))
               
  # test whether inside bounding rectangle
  xr <- Z$xrange
  yr <- Z$yrange
  eps <- sqrt(.Machine$double.eps)
  frameok <- (x >= xr[1] - eps) & (x <= xr[2] + eps) & 
             (y >= yr[1] - eps) & (y <= yr[2] + eps)
  
  if(!any(frameok)) {
    # all points OUTSIDE range - no further work needed
    if(!naok)
      warning("Internal error: all values NA")
    return(value)  # all NA
  }

  # consider only those points which are inside the frame
  xf <- x[frameok]
  yf <- y[frameok]
  # map locations to raster (row,col) coordinates
  if(strict)
    loc <- nearest.pixel(xf,yf,Z)
  else
    loc <- nearest.valid.pixel(xf,yf,Z)
  # look up image values
  vf <- Z$v[cbind(loc$row, loc$col)]
  
  # insert into answer
  value[frameok] <- vf

  if(!naok && any(is.na(value)))
    warning("Internal error: NA's generated")

  return(value)
}
  

rasterx.im <- function(x) {
  verifyclass(x, "im")
  v <- x$v
  xx <- x$xcol
  matrix(xx[col(v)], ncol=ncol(v), nrow=nrow(v))
}

rastery.im <- function(x) {
  verifyclass(x, "im")
  v <- x$v
  yy <- x$yrow
  matrix(yy[row(v)], ncol=ncol(v), nrow=nrow(v))
}

rasterxy.im <- function(x, drop=FALSE) {
  verifyclass(x, "im")
  v <- x$v
  xx <- x$xcol
  yy <- x$yrow
  if(!drop) {
    ans <- cbind(as.vector(xx[col(v)]), as.vector(yy[row(v)]))
  } else {
    ok <- !is.null(x$v)
    ans <- cbind(as.vector(xx[col(v)[ok]]), as.vector(yy[row(v)[ok]]))
  }
  colnames(ans) <- c("x", "y")
  return(ans)
}

##############

# methods for other functions

xtfrm.im <- function(x) { as.numeric(as.matrix.im(x)) }

as.matrix.im <- function(x, ...) {
  return(x$v)
}

as.array.im <- function(x, ...) {
  m <- as.matrix(x)
  a <- do.call(array, resolve.defaults(list(m),
                                       list(...),
                                       list(dim=c(dim(m), 1))))
  return(a)
}

as.data.frame.im <- function(x, ...) {
  verifyclass(x, "im")
  v <- x$v
  xx <- x$xcol[col(v)]
  yy <- x$yrow[row(v)]
  ok <- !is.na(v)
  xx <- as.vector(xx[ok])
  yy <- as.vector(yy[ok])
  # extract pixel values without losing factor info
  vv <- v[ok]
  dim(vv) <- NULL
  # 
  data.frame(x=xx, y=yy, value=vv, ...)
}
  
mean.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(mean(xvalues))
}

sum.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(sum(xvalues, ...))
}

median.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(median(xvalues, ...))
}

range.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(range(xvalues, ...))
}

max.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(max(xvalues, ...))
}

min.im <- function(x, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=TRUE]
  return(min(xvalues, ...))
}

hist.im <- function(x, ..., probability=FALSE) {
  xname <- short.deparse(substitute(x))
  verifyclass(x, "im")
  main <- paste("Histogram of", xname)
  # default plot arguments
  # extract pixel values
  values <- as.matrix(x)
  dim(values) <- NULL
  # barplot or histogram
  if(x$type %in% c("logical", "factor")) {
    # barplot
    tab <- table(values)
    probs <- tab/sum(tab)
    if(probability) {
      heights <- probs
      ylab <- "Probability"
    } else {
      heights <- tab
      ylab <- "Number of pixels"
    }
    mids <- do.call("barplot",
                   resolve.defaults(list(heights),
                                    list(...),
                                    list(xlab=paste("Pixel value"),
                                         ylab=ylab,
                                         main=main)))
    out <- list(counts=tab, probs=probs, heights=heights,
                mids=mids, xname=xname)
    class(out) <- "barplotdata"
  } else {
    # histogram
    values <- values[!is.na(values)]
    plotit <- resolve.defaults(list(...), list(plot=TRUE))$plot
    if(plotit) {
      ylab <- if(probability) "Probability density" else "Number of pixels"
      out <- do.call("hist.default",
                   resolve.defaults(list(values),
                                    list(...),
                                    list(probability=probability),
                                    list(xlab=paste("Pixel value"),
                                         ylab=ylab,
                                         main=main)))
      out$xname <- xname
    } else {
      # plot.default whinges if `probability' given when plot=FALSE
      out <- do.call("hist.default",
                   resolve.defaults(list(values),
                                    list(...)))
      # hack!
      out$xname <- xname
    }
  }
  return(invisible(out))
}

plot.barplotdata <- function(x, ...) {
  do.call("barplot",
          resolve.defaults(list(height=x$heights),
                           list(...),
                           list(main=paste("Histogram of ", x$xname))))
}

cut.im <- function(x, ...) {
  verifyclass(x, "im")
  vcut <- cut(as.numeric(as.matrix(x)), ...)
  return(im(vcut, xcol=x$xcol, yrow=x$yrow, unitname=unitname(x)))
}

quantile.im <- function(x, ...) {
  verifyclass(x, "im")
  q <- do.call("quantile",
               resolve.defaults(list(as.numeric(as.matrix(x))),
                                list(...),
                                list(na.rm=TRUE)))
  return(q)
}

integral.im <- function(x, ...) {
  verifyclass(x, "im")
  typ <- x$type
  if(!any(typ == c("integer", "real", "complex", "logical")))
    stop(paste("Don't know how to integrate an image of type", sQuote(typ)))
  a <- with(x, sum(v, na.rm=TRUE) * xstep * ystep)
  return(a)
}

conform.imagelist <- function(X, Zlist) {
  # determine points of X where all images in Zlist are defined
  ok <- rep.int(TRUE, length(X$x))
  for(i in seq_along(Zlist)) {
    Zi <- Zlist[[i]]
    ZiX <- Zi[X, drop=FALSE]
    ok <- ok & !is.na(ZiX)
  }
  return(ok)
}

split.im <- function(x, f, ..., drop=FALSE) {
  stopifnot(is.im(x))
  if(inherits(f, "tess")) 
    subsets <- tiles(f)
  else if(is.im(f)) {
    if(f$type != "factor")
      f <- eval.im(factor(f))
    subsets <- tiles(tess(image=f))
  } else stop("f should be a tessellation or a factor-valued image")
  if(!is.subset.owin(as.owin(x), as.owin(f)))
    stop("f does not cover the window of x")
  n <- length(subsets)
  out <- vector(mode="list", length=n)
  names(out) <- names(subsets)
  for(i in 1:n)
    out[[i]] <- x[subsets[[i]], drop=drop]
  if(drop)
    return(out)
  else 
    return(as.listof(out))
}

by.im <- function(data, INDICES, FUN, ...) {
  stopifnot(is.im(data))
  V <- split(data, INDICES)
  U <- lapply(V, FUN, ...)
  return(as.listof(U))
}

rebound.im <- function(x, rect) {
  stopifnot(is.im(x))
  stopifnot(is.owin(rect))
  rect <- as.rectangle(rect)
  stopifnot(is.subset.owin(as.rectangle(x), rect))
  # compute number of extra rows/columns
  dx <- x$xstep
  nleft  <- max(0, floor((x$xrange[1]-rect$xrange[1])/dx))
  nright <- max(0, floor((rect$xrange[2]-x$xrange[2])/dx))
  dy <- x$ystep
  nbot <- max(0, floor((x$yrange[1]-rect$yrange[1])/dy))
  ntop <- max(0, floor((rect$yrange[2]-x$yrange[2])/dy))
  # determine exact x and y ranges (to preserve original pixel locations)
  xrange.new <- x$xrange + c(-nleft, nright) * dx
  yrange.new <- x$yrange + c(-nbot,  ntop) * dy
  # expand pixel data matrix
  nr <- x$dim[1]
  nc <- x$dim[2]
  nrnew <- nbot  + nr + ntop
  ncnew <- nleft + nc + nright
  naval <- switch(x$type,
                  factor=,
                  integer=NA_integer_,
                  real=NA_real_,
                  character=NA_character_,
                  complex=NA_complex_,
                  NA)
  vnew <- matrix(naval, nrnew, ncnew)
  if(x$type != "factor") {
    vnew[nbot + (1:nr), nleft + (1:nc)] <- x$v
  } else {
    vnew[nbot + (1:nr), nleft + (1:nc)] <- as.integer(x$v)
    vnew <- factor(vnew, labels=levels(x))
    dim(vnew) <- c(nrnew, ncnew)
  }
  # build new image object
  xnew <- im(vnew,
             xrange = xrange.new,
             yrange = yrange.new,
             unitname = unitname(x))
  return(xnew)
}

sort.im <- function(x, ...) {
  verifyclass(x, "im")
  sort(as.vector(as.matrix(x)), ...)
}

dim.im <- function(x) { x$dim }

# colour images
rgbim <- function(R, G, B, maxColorValue=255) {
  eval.im(factor(rgbNA(as.vector(R), as.vector(G), as.vector(B),
                     maxColorValue=maxColorValue)))
}

hsvim <- function(H, S, V) {
  eval.im(factor(hsvNA(as.vector(H), as.vector(S), as.vector(V))))
}

scaletointerval <- function(x, from=0, to=1) {
  UseMethod("scaletointerval")
}

scaletointerval.default <- function(x, from=0, to=1) {
  rr <- range(x, na.rm=TRUE)
  b <- (to - from)/diff(rr)
  y <- from + b * (x - rr[1])
  return(y)
}

scaletointerval.im <- function(x, from=0, to=1) {
  v <- scaletointerval(x$v, from, to)
  y <- im(v, x$xcol, x$yrow, x$xrange, x$yrange, unitname(x))
  return(y)
}

zapsmall.im <- function(x, digits) {
  if(missing(digits))
    return(eval.im(zapsmall(x)))
  return(eval.im(zapsmall(x, digits=digits)))
}

