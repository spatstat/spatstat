#
#  quadratcount.R
#
#  $Revision: 1.54 $  $Date: 2016/02/16 01:39:12 $
#

quadratcount <- function(X, ...) {
  UseMethod("quadratcount")
}

quadratcount.splitppp <- function(X, ...) {
  solapply(X, quadratcount, ...)
}

quadratcount.ppp <- function(X, nx=5, ny=nx, ...,
                             xbreaks=NULL, ybreaks=NULL,
                             tess=NULL)  {
  verifyclass(X, "ppp")
  W <- X$window

  if(is.null(tess)) {
    # rectangular boundaries 
    if(!is.numeric(nx))
      stop("nx should be numeric")
    # start with rectangular tessellation
    tess <- quadrats(as.rectangle(W),
                     nx=nx, ny=ny, xbreaks=xbreaks, ybreaks=ybreaks)
    # fast code for counting points in rectangular grid
    Xcount <- rectquadrat.countEngine(X$x, X$y, tess$xgrid, tess$ygrid)
    #
    if(W$type != "rectangle") {
      # intersections of rectangles with window including empty intersections
      tess <- quadrats(X,
                       nx=nx, ny=ny, xbreaks=xbreaks, ybreaks=ybreaks,
                       keepempty=TRUE)
      # now delete the empty quadrats and the corresponding counts
      nonempty <- !tiles.empty(tess)
#     WAS: nonempty <- !unlist(lapply(tiles(tess), is.empty))
      if(!any(nonempty))
        stop("All tiles are empty")
      if(!all(nonempty)) {
#        ntiles <- sum(nonempty)
        tess   <- tess[nonempty]
        Xcount <- t(Xcount)[nonempty]
        # matrices and tables are in row-major order,
        # tiles in a rectangular tessellation are in column-major order
        Xcount <- array(Xcount,
                        dimnames=list(tile=tilenames(tess)))
        class(Xcount) <- "table"
      }
    }
  } else {
    # user-supplied tessellation
    if(!inherits(tess, "tess")) {
      tess <- try(as.tess(tess), silent=TRUE)
      if(inherits(tess, "try-error"))
        stop("The argument tess should be a tessellation", call.=FALSE)
    }
    if(tess$type == "rect") {
      # fast code for counting points in rectangular grid
      Xcount <- rectquadrat.countEngine(X$x, X$y, tess$xgrid, tess$ygrid)
    } else {
      # quadrats are another type of tessellation
      Y <- cut(X, tess)
      if(anyNA(marks(Y)))
        warning("Tessellation does not contain all the points of X")
      Xcount <- table(tile=marks(Y))
    }
  }
  attr(Xcount, "tess") <- tess
  class(Xcount) <- c("quadratcount", class(Xcount))
  return(Xcount)
}

plot.quadratcount <- function(x, ...,
                              add=FALSE, entries=as.vector(t(as.table(x))),
                              dx=0, dy=0, show.tiles=TRUE,
                              textargs = list()) {
  xname <- short.deparse(substitute(x))
  tess <- attr(x, "tess")
  # add=FALSE, show.tiles=TRUE  => plot tiles + numbers
  # add=FALSE, show.tiles=FALSE => plot window (add=FALSE) + numbers
  # add=TRUE,  show.tiles=TRUE  => plot tiles  (add=TRUE) + numbers
  # add=TRUE,  show.tiles=FALSE => plot numbers
  if(show.tiles || !add) {
    context <- if(show.tiles) tess else as.owin(tess)
    do.call(plot,
            resolve.defaults(list(context, add=add),
                             list(...),
                             list(main=xname),
                             .StripNull=TRUE))
  }
  if(!is.null(entries)) {
    labels <- paste(as.vector(entries))
    til <- tiles(tess)
    incircles <- lapply(til, incircle)
    x0 <- sapply(incircles, getElement, name="x")
    y0 <- sapply(incircles, getElement, name="y")
    ra <- sapply(incircles, getElement, name="r")
    do.call.matched(text.default,
                    resolve.defaults(list(x=x0 + dx * ra, y = y0 + dy * ra),
                                     list(labels=labels),
                                     textargs, 
                                     list(...)))
  }
  return(invisible(NULL))
}

rectquadrat.breaks <- function(xr, yr, nx=5, ny=nx, xbreaks=NULL, ybreaks=NULL) {
  if(is.null(xbreaks))
    xbreaks <- seq(from=xr[1], to=xr[2], length.out=nx+1)
  else if(min(xbreaks) > xr[1] || max(xbreaks) < xr[2])
    stop("xbreaks do not span the range of x coordinates in the window")
  if(is.null(ybreaks))
    ybreaks <- seq(from=yr[1], to=yr[2], length.out=ny+1)
  else if(min(ybreaks) > yr[1] || max(ybreaks) < yr[2])
    stop("ybreaks do not span the range of y coordinates in the window")
  return(list(xbreaks=xbreaks, ybreaks=ybreaks))
}

rectquadrat.countEngine <- function(x, y, xbreaks, ybreaks, weights) {
  if(length(x) > 0) {
    # check validity of breaks
    if(!all(inside.range(range(x), range(xbreaks))))
      stop("xbreaks do not span the actual range of x coordinates in data")
    if(!all(inside.range(range(y), range(ybreaks))))
      stop("ybreaks do not span the actual range of y coordinates in data")
  }
  # WAS: 
  # xg <- cut(x, breaks=xbreaks, include.lowest=TRUE)
  # yg <- cut(y, breaks=ybreaks, include.lowest=TRUE)
  xg <- fastFindInterval(x, xbreaks, labels=TRUE)
  yg <- fastFindInterval(y, ybreaks, labels=TRUE)
  if(missing(weights)) 
    sumz <- table(list(y=yg, x=xg))
  else {
    sumz <- tapply(weights, list(y=yg, x=xg), sum)
    if(any(nbg <- is.na(sumz)))
      sumz[nbg] <- 0
  }
  # reverse order of y 
  sumz <- sumz[rev(seq_len(nrow(sumz))), ]
  sumz <- as.table(sumz)
  #
  attr(sumz, "xbreaks") <- xbreaks
  attr(sumz, "ybreaks") <- ybreaks
  return(sumz)
}

quadrats <- function(X, nx=5, ny=nx, xbreaks = NULL, ybreaks = NULL,
                     keepempty=FALSE) {
  W <- as.owin(X)
  xr <- W$xrange
  yr <- W$yrange
  b <- rectquadrat.breaks(xr, yr, nx, ny, xbreaks, ybreaks)
  # rectangular tiles
  Z <- tess(xgrid=b$xbreaks, ygrid=b$ybreaks, unitname=unitname(W))
  if(W$type != "rectangle") {
    # intersect rectangular tiles with window W
    if(!keepempty) {
      Z <- intersect.tess(Z, W)
    } else {
      til <- tiles(Z)
      for(i in seq_along(til))
        til[[i]] <- intersect.owin(til[[i]], W)
      Z <- tess(tiles=til, window=W, keepempty=TRUE)
    }
  }
  return(Z)
}

as.tess.quadratcount <- function(X) {
  return(attr(X, "tess"))
}

as.owin.quadratcount <- function(W, ..., fatal=TRUE) {
  return(as.owin(as.tess(W), ..., fatal=fatal))
}

domain.quadratcount <- Window.quadratcount <- function(X, ...) { as.owin(X) }

intensity.quadratcount <- function(X, ..., image=FALSE) {
  Y <- as.tess(X)
  a <- tile.areas(Y)
  ## in the rectangular case, tiles are indexed in column-major order
  if(Y$type == "rect" && length(dim(X)) > 1) 
    a <- matrix(a, byrow=TRUE, nrow(X), ncol(X))
  lambda <- X/a
  if(!image) {
    trap.extra.arguments(...)
    class(lambda) <- "table"
    attr(lambda, "tess") <- NULL
    return(lambda)
  }
  ## again to handle rectangular case
  lambda <- as.vector(t(lambda))
  tileid <- as.im(Y, ...)
  result <- eval.im(lambda[tileid])
  return(result)
}

## The shift method is undocumented.
## It is only needed in plot.listof / plot.solist / plot.layered

shift.quadratcount <- function(X, ...) {
  attr(X, "tess") <- te <- shift(attr(X, "tess"), ...)
  attr(X, "lastshift") <- getlastshift(te)
  return(X)
}

