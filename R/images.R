#
#       images.R
#
#      $Revision: 1.147 $     $Date: 2018/04/12 10:05:05 $
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
  if(!is.null(dim(mat))) {
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
    xcol <- seq(from=xrange[1L] + xstep/2, to=xrange[2L] - xstep/2, length.out=nc)
  } else if(length(xcol) > 1) {
    # use 'xcol'
    # ensure spacing is constant
    xcol <- seq(from=min(xcol), to=max(xcol), length.out=length(xcol))
    xstep <- diff(xcol)[1L]
    xrange <- range(xcol) + c(-1,1) * xstep/2
  } else stop("Cannot determine pixel width")
  
  if((miss.yrow || length(yrow) <= 1) && !is.null(yrange)) {
    # use 'yrange'
    ystep <- diff(yrange)/nr
    yrow <- seq(from=yrange[1L] + ystep/2, to=yrange[2L] - ystep/2, length.out=nr)
  } else if(length(yrow) > 1) {
    # use 'yrow'
    # ensure spacing is constant
    yrow <- seq(from=min(yrow), to=max(yrow), length.out=length(yrow))
    ystep <- diff(yrow)[1L]
    yrange <- range(yrow) + c(-1,1) * ystep/2
  }  else stop("Cannot determine pixel height")

  unitname <- as.unitname(unitname)

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
                   bottomleft={ c(W$xrange[1L], W$yrange[1L]) })
    return(shift(X, -locn))
  }
  X$xrange <- X$xrange + vec[1L]
  X$yrange <- X$yrange + vec[2L]
  X$xcol <- X$xcol + vec[1L]
  X$yrow <- X$yrow + vec[2L]
  attr(X, "lastshift") <- vec
  return(X)
}

"Frame<-.im" <- function(X, value) {
  stopifnot(is.rectangle(value))
  if(!is.subset.owin(value, Frame(X))) {
    ## first expand
    X <- X[value, drop=FALSE]
  }
  X[value, drop=TRUE]
}


"[.im" <- local({

  disjoint <- function(r, s) { (r[2L] < s[1L]) || (r[1L] > s[2L])  }
  clip <- function(r, s) { c(max(r[1L],s[1L]), min(r[2L],s[2L])) }
  inrange <- function(x, r) { (x >= r[1L]) & (x <= r[2L]) }

  Extract.im <- function(x, i, j, ...,
                         drop=TRUE, tight=FALSE, raster=NULL,
                         rescue=is.owin(i)) {

    ## detect 'blank' arguments like second argument in x[i, ] 
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
      ## no indices: return entire image 
      out <- if(is.null(raster)) x else as.im(raster)
      xy <- expand.grid(y=out$yrow,x=out$xcol)
      if(!is.null(raster)) {
        ## resample image on new pixel raster
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
      ## .................................................................
      ## Try spatial index
      ## .................................................................
      if(verifyclass(i, "owin", fatal=FALSE)) {

        if(jtype == "given")
          warning("Argument j ignored")
      
        ## 'i' is a window
        ## if drop = FALSE, just set values outside window to NA
        ## if drop = TRUE, extract values for all pixels inside window
        ##                 as an image (if 'i' is a rectangle)
        ##                 or as a vector (otherwise)

        ## determine pixel raster for output
        if(!is.null(raster)) {
          out <- as.im(raster)
          do.resample <- TRUE
        } else if(is.subset.owin(i, as.owin(x))) {
          out <- x
          do.resample <- FALSE
        } else {
          ## new window does not contain data window: expand it
          bb <- boundingbox(as.rectangle(i), as.rectangle(x))
          rr <- if(is.mask(i)) i else x
          xcol <- prolongseq(rr$xcol, bb$xrange, rr$xstep)
          yrow <- prolongseq(rr$yrow, bb$yrange, rr$ystep)
          out <- list(xcol=xcol, yrow=yrow)
          do.resample <- TRUE
        }
        xy <- expand.grid(y=out$yrow,x=out$xcol)
        if(do.resample) {
          ## resample image on new pixel raster
          values <- lookup.im(x, xy$x, xy$y, naok=TRUE)
          out <- im(values, out$xcol, out$yrow, unitname=unitname(out))
        }
        inside <- inside.owin(xy$x, xy$y, i)
        if(!drop) {
          ## set other pixels to NA and return image
          out$v[!inside] <- NA
          if(!tight)
            return(out)
        } else if(!(rescue && i$type == "rectangle")) {
          ## return pixel values
          values <- out$v[inside]
          return(values)
        }
        ## return image in smaller rectangle
        if(disjoint(i$xrange, x$xrange) || disjoint(i$yrange, x$yrange))
          ## empty intersection
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
        result <- do.call(im, c(marg, xarg, yarg))
        return(result)
      }
      if(verifyclass(i, "im", fatal=FALSE)) {
        if(jtype == "given")
          warning("Argument j ignored")
        ## logical images OK
        if(i$type == "logical") {
          ## convert to window
          w <- as.owin(eval.im(ifelse1NA(i)))
          return(x[w, drop=drop, ..., raster=raster])
        } else stop("Subset argument \'i\' is an image, but not of logical type")
      }

      if(inherits(i, "linnet")) {
        #' linear network
        if(jtype == "given")
          warning("Argument j ignored")
        W <- raster %orifnull% as.owin(x)
        M <- as.mask.psp(as.psp(i), W=W, ...)
        xM <- x[M, drop=drop]
        if(is.im(xM)) xM <- linim(i, xM)
        return(xM)
      }
      
      if(is.ppp(i)) {
        ## 'i' is a point pattern 
        if(jtype == "given")
          warning("Argument j ignored")
        ## Look up the greyscale values for the points of the pattern
        values <- lookup.im(x, i$x, i$y, naok=TRUE)
        if(drop) 
          values <- values[!is.na(values)]
        if(length(values) == 0) 
          ## ensure the zero-length vector is of the right type
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
    ## ............... not a spatial index .............................

    ## Try indexing as a matrix

    ## Construct a matrix index call for possible re-use
    M <- as.matrix(x)
    ## suppress warnings from code checkers
    dont.complain.about(M)
    ##
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
    ## try it
    y <- try(eval(as.call(ycall)), silent=TRUE)
    if(!inherits(y, "try-error")) {
      ## valid subset index for a matrix
      if(rescue) {
        ## check whether it's a rectangular block, in correct order
        RR <- row(x$v)
        CC <- col(x$v)
        rcall <- ycall
        rcall[[2L]] <- quote(RR)
        ccall <- ycall
        ccall[[2L]] <- quote(CC)
        rr <- eval(as.call(rcall))
        cc <- eval(as.call(ccall))
        rseq <- sort(unique(as.vector(rr)))
        cseq <- sort(unique(as.vector(cc)))
        if(all(diff(rseq) == 1) && all(diff(cseq) == 1) &&
           (length(rr) == length(rseq) * length(cseq)) &&
           all(rr == RR[rseq, cseq]) && all(cc == CC[rseq,cseq])) {
          ## yes - make image
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
      ## return pixel values (possibly as matrix)
      return(y)
    }

    ## Last chance!
    if(itype == "given" &&
       !is.matrix(i) &&
       !is.null(ip <- as.ppp(i, W=as.owin(x), fatal=FALSE, check=FALSE))) {
      ## 'i' is convertible to a point pattern 
      ## Look up the greyscale values for the points of the pattern
      values <- lookup.im(x, ip$x, ip$y, naok=TRUE)
      if(drop) 
        values <- values[!is.na(values)]
      if(length(values) == 0) 
        ## ensure the zero-length vector is of the right type
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

  Extract.im
})

update.im <- function(object, ...) {
  ## update internal structure of image after manipulation
  X <- object
  mat <- X$v
  typ <- typeof(mat)
  if(typ == "double")
    typ <- "real"
  ## deal with factor case
  if(is.factor(mat)) {
    typ <- "factor"
  } else if(!is.null(lev <- levels(mat))) {
    typ <- "factor"
    X$v <- factor(mat, levels=lev)
  }
  X$type <- typ
  return(X)
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

  stopifnot(is.im(value) || is.vector(value) ||
            is.matrix(value) || is.array(value) || is.factor(value))
  if(is.im(value)) 
    value <- value$v

  if(itype == "missing" && jtype == "missing") {
    # no index provided
    # set all pixels to 'value'
    v <- X$v
    if(!is.factor(value)) {
      v[!is.na(v)] <- value
    } else {
      vnew <- matrix(NA_integer_, ncol(v), nrow(v))
      vnew[!is.na(v)] <- as.integer(value)
      v <- factor(vnew, labels=levels(value))
    }
    X$v <- v
    return(update(X))
  }
  if(itype == "given") {
    # ..................... Try a spatial index ....................
    if(verifyclass(i, "owin", fatal=FALSE)) {
      if(jtype == "given") warning("Index j ignored")
      # 'i' is a window
      if(is.empty(i))
        return(X)
      rxy <- rasterxy.mask(W)
      xx <- rxy$x
      yy <- rxy$y
      ok <- inside.owin(xx, yy, i)
      X$v[ok] <- value
      X$type <- ifelse(is.factor(X$v), "factor", typeof(X$v))
      return(update(X))
    }
    if(verifyclass(i, "im", fatal=FALSE) && i$type == "logical") {
      if(jtype == "given") warning("Index j ignored")
      # convert logical vector to window where entries are TRUE
      i <- as.owin(eval.im(ifelse1NA(i)))
      # continue as above
      rxy <- rasterxy.mask(W)
      xx <- rxy$x
      yy <- rxy$y
      ok <- inside.owin(xx, yy, i)
      X$v[ok] <- value
      X$type <- ifelse(is.factor(X$v), "factor", typeof(X$v))
      return(update(X))
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
      X$type <- ifelse(is.factor(X$v), "factor", typeof(X$v))
      return(update(X))
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
  if(!inherits(litmus, "try-error")){
    X$type <- ifelse(is.factor(X$v), "factor", typeof(X$v))
    return(update(X))
  }
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
    X$type <- ifelse(is.factor(X$v), "factor", typeof(X$v))
    return(update(X))
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

nearest.pixel <- function(x,y, Z) {
  stopifnot(is.im(Z) || is.mask(Z))
  if(length(x) > 0) {
    nr <- Z$dim[1L]
    nc <- Z$dim[2L]
    cc <- round(1 + (x - Z$xcol[1L])/Z$xstep)
    rr <- round(1 + (y - Z$yrow[1L])/Z$ystep)
    cc <- pmax.int(1,pmin.int(cc, nc))
    rr <- pmax.int(1,pmin.int(rr, nr))
  } else cc <- rr <- integer(0)
  return(list(row=rr, col=cc))
}

# Explores the 3 x 3 neighbourhood of nearest.pixel
# and finds the nearest pixel that is not NA

nearest.valid.pixel <- function(x, y, Z) {
  rc <- nearest.pixel(x,y,Z) # checks that Z is an 'im' or 'mask'
  rr <- rc$row
  cc <- rc$col
  # check whether any pixels are outside image domain
  inside <- as.owin(Z)$m
  miss <- !inside[cbind(rr, cc)]
  if(!any(miss))
    return(rc)
  # for offending pixels, explore 3 x 3 neighbourhood
  nr <- Z$dim[1L]
  nc <- Z$dim[2L]
  xcol <- Z$xcol
  yrow <- Z$yrow
  for(i in which(miss)) {
    rows <- rr[i] + c(-1L,0L,1L)
    cols <- cc[i] + c(-1L,0L,1L)
    rows <- unique(pmax.int(1, pmin.int(rows, nr)))
    cols <- unique(pmax.int(1, pmin.int(cols, nc)))
    rcp <- expand.grid(row=rows, col=cols)
    ok <- inside[as.matrix(rcp)]
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
  frameok <- (x >= xr[1L] - eps) & (x <= xr[2L] + eps) & 
             (y >= yr[1L] - eps) & (y <= yr[2L] + eps)
  
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

  if(!naok && anyNA(value))
    warning("Internal error: NA's generated")

  return(value)
}
  

## low level

rasterx.im <- function(x) {
  verifyclass(x, "im")
  xx <- x$xcol
  matrix(xx[col(x)], ncol=ncol(x), nrow=nrow(x))
}

rastery.im <- function(x) {
  verifyclass(x, "im")
  yy <- x$yrow
  matrix(yy[row(x)], ncol=ncol(x), nrow=nrow(x))
}

rasterxy.im <- function(x, drop=FALSE) {
  verifyclass(x, "im")
  xx <- x$xcol
  yy <- x$yrow
  ans <- cbind(x=as.vector(xx[col(x)]),
               y=as.vector(yy[row(x)]))
  if(drop) {
    ok <- as.vector(!is.na(x$v))
    ans <- ans[ok, , drop=FALSE]
  }
  return(ans)
}

## user interface 

raster.x <- function(w, drop=FALSE) {
  if(is.owin(w)) return(rasterx.mask(w, drop=drop))
  if(!is.im(w)) stop("w should be a window or an image")
  x <- w$xcol[col(w)]
  x <- if(drop) x[!is.na(w$v), drop=TRUE] else array(x, dim=w$dim)
  return(x)
}
  
raster.y <- function(w, drop=FALSE) {
  if(is.owin(w)) return(rastery.mask(w, drop=drop))
  if(!is.im(w)) stop("w should be a window or an image")
  y <- w$yrow[row(w)]
  y <- if(drop) y[!is.na(w$v), drop=TRUE] else array(y, dim=w$dim)
  return(y)
}

raster.xy <- function(w, drop=FALSE) {
  if(is.owin(w)) return(rasterxy.mask(w, drop=drop))
  if(!is.im(w)) stop("w should be a window or an image")
  x <- w$xcol[col(w)]
  y <- w$yrow[row(w)]
  if(drop) {
    ok <- !is.na(w$v)
    x <- x[ok, drop=TRUE]
    y <- y[ok, drop=TRUE]
  }
  return(list(x=as.numeric(x),
              y=as.numeric(y)))
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

mean.im <- function(x, trim=0, na.rm=TRUE, ...) {
  verifyclass(x, "im")
  xvalues <- x[drop=na.rm]
  return(mean(xvalues, trim=trim, na.rm=na.rm))
}

## arguments of generic 'median' will change in R 3.4
median.im <- if("..." %in% names(formals(median))) {
   function(x, na.rm=TRUE, ...) {
    verifyclass(x, "im")
    xvalues <- x[drop=na.rm]
    return(median(xvalues, ...))
  }
} else {
   function(x, na.rm=TRUE) {
    verifyclass(x, "im")
    xvalues <- x[drop=na.rm]
    return(median(xvalues))
  }
}

where.max <- function(x, first=TRUE) {
  stopifnot(is.im(x))
  if(first) { 
    ## find the first maximum
    v <- x$v
    locn <- which.max(as.vector(v))  # ignores NA, NaN
    locrow <- as.vector(row(v))[locn]
    loccol <- as.vector(col(v))[locn]
  } else {
    ## find all maxima
    xmax <- max(x)
    M <- solutionset(x == xmax)
    loc <- which(M$m, arr.ind=TRUE)
    locrow <- loc[,1L]
    loccol <- loc[,2L]
  }
  xx <- x$xcol[loccol]
  yy <- x$yrow[locrow]
  return(ppp(x=xx, y=yy, window=Window(x)))
}

where.min <- function(x, first=TRUE) {
  stopifnot(is.im(x))
  if(first) { 
    ## find the first minimum
    v <- x$v
    locn <- which.min(as.vector(v))  # ignores NA, NaN
    locrow <- as.vector(row(v))[locn]
    loccol <- as.vector(col(v))[locn]
  } else {
    ## find all minima
    xmin <- min(x)
    M <- solutionset(x == xmin)
    loc <- which(M$m, arr.ind=TRUE)
    locrow <- loc[,1L]
    loccol <- loc[,2L]
  }
  xx <- x$xcol[loccol]
  yy <- x$yrow[locrow]
  return(ppp(x=xx, y=yy, window=Window(x)))
}

## the following ensures that 'sd' works

as.double.im <- function(x, ...) { as.double(x[], ...) }

##

hist.im <- function(x, ..., probability=FALSE, xname) {
  if(missing(xname) || is.null(xname)) xname <- short.deparse(substitute(x))
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
    mids <- do.call(barplot,
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
      out <- do.call(hist.default,
                     resolve.defaults(list(values),
                                      list(...),
                                      list(freq=!probability,
                                           xlab="Pixel value",
                                           ylab=ylab,
                                           main=main)))
      out$xname <- xname
    } else {
      # plot.default whinges if `probability' given when plot=FALSE
      out <- do.call(hist.default,
                   resolve.defaults(list(values),
                                    list(...)))
      # hack!
      out$xname <- xname
    }
  }
  return(invisible(out))
}

plot.barplotdata <- function(x, ...) {
  do.call(barplot,
          resolve.defaults(list(height=x$heights),
                           list(...),
                           list(main=paste("Histogram of ", x$xname))))
}

cut.im <- function(x, ...) {
  verifyclass(x, "im")
  typ <- x$type
  if(typ %in% c("factor", "logical", "character")) 
    stop(paste0("cut.im is not defined for ", typ, "-valued images"),
         call.=FALSE)
  vcut <- cut(as.numeric(as.matrix(x)), ...)
  return(im(vcut,
            xcol=x$xcol, yrow=x$yrow,
            xrange=x$xrange, yrange=x$yrange,
            unitname=unitname(x)))
}

quantile.im <- function(x, ...) {
  verifyclass(x, "im")
  q <- do.call(quantile,
               resolve.defaults(list(as.numeric(as.matrix(x))),
                                list(...),
                                list(na.rm=TRUE)))
  return(q)
}

integral <- function(f, domain=NULL, ...) {
  UseMethod("integral")
}

integral.im <- function(f, domain=NULL, ...) {
  verifyclass(f, "im")
  typ <- f$type
  if(!any(typ == c("integer", "real", "complex", "logical")))
    stop(paste("Don't know how to integrate an image of type", sQuote(typ)))
  if(!is.null(domain)) {
    if(is.tess(domain)) return(sapply(tiles(domain), integral.im, f=f))
    f <- f[domain, drop=FALSE, tight=TRUE]
  }
  a <- with(f, sum(v, na.rm=TRUE) * xstep * ystep)
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
    return(as.solist(out))
}

by.im <- function(data, INDICES, FUN, ...) {
  stopifnot(is.im(data))
  V <- split(data, INDICES)
  U <- lapply(V, FUN, ...)
  return(as.solist(U, demote=TRUE))
}

rebound.im <- function(x, rect) {
  stopifnot(is.im(x))
  stopifnot(is.owin(rect))
  rect <- as.rectangle(rect)
  stopifnot(is.subset.owin(as.rectangle(x), rect))
  # compute number of extra rows/columns
  dx <- x$xstep
  nleft  <- max(0, floor((x$xrange[1L]-rect$xrange[1L])/dx))
  nright <- max(0, floor((rect$xrange[2L]-x$xrange[2L])/dx))
  dy <- x$ystep
  nbot <- max(0, floor((x$yrange[1L]-rect$yrange[1L])/dy))
  ntop <- max(0, floor((rect$yrange[2L]-x$yrange[2L])/dy))
  # determine exact x and y ranges (to preserve original pixel locations)
  xrange.new <- x$xrange + c(-nleft, nright) * dx
  yrange.new <- x$yrange + c(-nbot,  ntop) * dy
  # expand pixel data matrix
  nr <- x$dim[1L]
  nc <- x$dim[2L]
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
rgbim <- function(R, G, B, A=NULL, maxColorValue=255, autoscale=FALSE) {
  if(autoscale) {
    R <- scaletointerval(R, 0, maxColorValue)
    G <- scaletointerval(G, 0, maxColorValue)
    B <- scaletointerval(B, 0, maxColorValue)
    if(!is.null(A))
      A <- scaletointerval(A, 0, maxColorValue)
  }
  Z <- eval.im(factor(rgbNA(as.vector(R), as.vector(G), as.vector(B),
                            as.vector(A),
                            maxColorValue=maxColorValue)))
  return(Z)
}

hsvim <- function(H, S, V, A=NULL, autoscale=FALSE) {
  if(autoscale) {
    H <- scaletointerval(H, 0, 1)
    S <- scaletointerval(S, 0, 1)
    V <- scaletointerval(V, 0, 1)
    if(!is.null(A))
      A <- scaletointerval(A, 0, 1)
  }
  Z <- eval.im(factor(hsvNA(as.vector(H), as.vector(S), as.vector(V),
                            as.vector(A))))
  return(Z)
}

scaletointerval <- function(x, from=0, to=1, xrange=range(x)) {
  UseMethod("scaletointerval")
}

scaletointerval.default <- function(x, from=0, to=1, xrange=range(x)) {
  x <- as.numeric(x)
  rr <- if(missing(xrange)) range(x, na.rm=TRUE) else as.numeric(xrange)
  b <- as.numeric(to - from)/diff(rr)
  if(is.finite(b)) {
    y <- from + b * (x - rr[1L])
  } else {
    y <- (from+to)/2 + 0 * x
  }
  y[] <- pmin(pmax(y[], from), to)
  return(y)
}

scaletointerval.im <- function(x, from=0, to=1, xrange=range(x)) {
  v <- scaletointerval(x$v, from, to, xrange=xrange)
  y <- im(v, x$xcol, x$yrow, x$xrange, x$yrange, unitname(x))
  return(y)
}

zapsmall.im <- function(x, digits) {
  if(missing(digits))
    return(eval.im(zapsmall(x)))
  return(eval.im(zapsmall(x, digits=digits)))
}

domain.im <- Window.im <- function(X, ...) { as.owin(X) }

"Window<-.im" <- function(X, ..., value) {
  verifyclass(value, "owin")
  X[value, drop=FALSE]
}

padimage <- function(X, value=NA, n=1, W=NULL) {
  stopifnot(is.im(X))
  stopifnot(length(value) == 1)
  if(!missing(n) && !is.null(W)) 
    stop("Arguments n and W are incompatible", call.=FALSE)
  padW <- !is.null(W)
  if(isfac <- (X$type == "factor")) {
    ## handle factors
    levX <- levels(X)
    if(is.factor(value)) {
      stopifnot(identical(levels(X), levels(value)))
    } else {
      value <- factor(value, levels=levX)
    }
    X <- eval.im(as.integer(X))
    value <- as.integer(value)
  }
  if(!padW) {
    ## pad by 'n' pixels
    nn <- rep(n, 4)
    nleft   <- nn[1L]
    nright  <- nn[2L]
    nbottom <- nn[3L]
    ntop    <- nn[4L]
  } else {
    ## pad out to window W
    FX <- Frame(X)
    B <- boundingbox(Frame(W), FX)
    nleft   <- max(1, round((FX$xrange[1L] - B$xrange[1L])/X$xstep))
    nright  <- max(1, round((B$xrange[2L] - FX$xrange[2L])/X$xstep))
    nbottom <- max(1, round((FX$yrange[1L] - B$yrange[1L])/X$ystep))
    ntop    <- max(1, round((B$yrange[2L] - FX$yrange[2L])/X$ystep))
  }
  mX <- as.matrix(X)
  dd <- dim(mX)
  mX <- cbind(matrix(value, dd[1L], nleft, byrow=TRUE),
              as.matrix(X),
              matrix(value, dd[1L], nright, byrow=TRUE))
  dd <- dim(mX)
  mX <- rbind(matrix(rev(value), nbottom, dd[2L]),
              mX,
              matrix(value, ntop, dd[2L]))
  xcol <- with(X,
               c(xcol[1L]     - (nleft:1) * xstep,
                 xcol,
                 xcol[length(xcol)] + (1:nright) * xstep))
  yrow <- with(X,
               c(yrow[1L]     - (nbottom:1) * ystep,
                 yrow,
                 yrow[length(yrow)] + (1:ntop) * ystep))
  xr <- with(X, xrange + c(-nleft, nright) * xstep)
  yr <- with(X, yrange + c(-nbottom, ntop) * ystep)
  Y <- im(mX,
          xcol=xcol, yrow=yrow, xrange=xr, yrange=yr,
          unitname=unitname(X))
  if(isfac)
    Y <- eval.im(factor(Y, levels=seq_along(levX), labels=levX))
  if(padW && !is.rectangle(W)) 
    Y <- Y[W, drop=FALSE]
  return(Y)
}

as.function.im <- function(x, ...) {
  Z <- x
  f <- function(x,y) { Z[list(x=x, y=y)] }
  g <- funxy(f, Window(x))
  return(g)
}

anyNA.im <- function(x, recursive=FALSE) {
  anyNA(x$v)
}

