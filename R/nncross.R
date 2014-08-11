#
#   nncross.R
#
#
#    $Revision: 1.18 $  $Date: 2012/10/13 02:02:00 $
#
#  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2012
#  Licence: GNU Public Licence >= 2

nncross <- function(X, Y, ...) {
  UseMethod("nncross")
}

nncross.default <- function(X, Y, ...) {
  X <- as.ppp(X, W=bounding.box.xy)
  nncross(X, Y, ...)
}

nncross.ppp <- function(X, Y, iX=NULL, iY=NULL,
                    what = c("dist", "which"),
                    ...,
                    sortby=c("range", "var", "x", "y"),
                    is.sorted.X = FALSE,
                    is.sorted.Y = FALSE) {
  stopifnot(is.ppp(Y) || is.psp(Y))
  sortby <- match.arg(sortby)
  what   <- match.arg(what, choices=c("dist", "which"), several.ok=TRUE)

  # deal with null cases
  nX <- X$n
  nY <- Y$n
  if(nX == 0)
    return(data.frame(dist=numeric(0), which=integer(0))[, what])
  if(nY == 0)
    return(data.frame(dist=rep(Inf, nX), which=rep(NA, nX))[, what])

  # Y is a line segment pattern 
  if(is.psp(Y))
    return(ppllengine(X,Y,"distance")[, what])

  # Y is a point pattern
  if(is.null(iX) != is.null(iY))
    stop("If one of iX, iY is given, then both must be given")
  exclude <- (!is.null(iX) || !is.null(iY))
  if(exclude) {
    stopifnot(is.integer(iX) && is.integer(iY))
    if(length(iX) != nX)
      stop("length of iX does not match the number of points in X")
    if(length(iY) != nY)
      stop("length of iY does not match the number of points in Y")
  }

  if((is.sorted.X || is.sorted.Y) && !(sortby %in% c("x", "y")))
     stop(paste("If data are already sorted,",
                "the sorting coordinate must be specified explicitly",
                "using sortby = \"x\" or \"y\""))

  # decide whether to sort on x or y coordinate
  switch(sortby,
         range = {
           WY <- as.owin(Y)
           sortby.y <- (diff(WY$xrange) < diff(WY$yrange))
         },
         var = {
           sortby.y <- (var(Y$x) < var(Y$y))
         },
         x={ sortby.y <- FALSE},
         y={ sortby.y <- TRUE}
         )

  # The C code expects points to be sorted by y coordinate.
  if(sortby.y) {
    Xx <- X$x
    Xy <- X$y
    Yx <- Y$x
    Yy <- Y$y
  } else {
    Xx <- X$y
    Xy <- X$x
    Yx <- Y$y
    Yy <- Y$x
  }
  # sort only if needed
  if(!is.sorted.X){
    oX <- fave.order(Xy)
    Xx <- Xx[oX]
    Xy <- Xy[oX]
    if(exclude) iX <- iX[oX]
  }
  if (!is.sorted.Y){
    oY <- fave.order(Yy)
    Yx <- Yx[oY]
    Yy <- Yy[oY]
    if(exclude) iY <- iY[oY]
  }
  
  # call C code
  want.dist  <- "dist" %in% what 
  want.which <- "which" %in% what
  want.both  <- want.dist && want.which
  Cfun <- paste("nnX",
                if(exclude) "E" else "",
                if(want.both) "" else if(want.dist) "dist" else "which",
                sep="")
  nndv <- if(want.dist) numeric(nX) else numeric(1)
  nnwh <- if(want.which) integer(nX) else integer(1)
  if(!exclude) iX <- iY <- integer(1)
  DUP <- spatstat.options("dupC")
  huge <- 1.1 * diameter(bounding.box(as.rectangle(X), as.rectangle(Y)))
  
  z <- .C(Cfun,
          n1=as.integer(nX),
          x1=as.double(Xx),
          y1=as.double(Xy),
          id1=as.integer(iX),
          n2=as.integer(nY),
          x2=as.double(Yx),
          y2=as.double(Yy),
          id2=as.integer(iY),
          nnd=as.double(nndv),
          nnwhich=as.integer(nnwh),
          huge=as.double(huge),
          DUP=DUP,
          PACKAGE="spatstat")

  if(want.which) {
    nnwcode <- z$nnwhich #sic. C code now increments by 1
    if(any(uhoh <- (nnwcode == 0))) {
      warning("NA's produced in nncross()$which")
      nnwcode[uhoh] <- NA
    }
  }
  
  # reinterpret in original ordering
  if(is.sorted.X){
    if(want.dist) nndv <- z$nnd
    if(want.which) nnwh <- if(is.sorted.Y) nnwcode else oY[nnwcode]
  } else {
    if(want.dist) nndv[oX] <- z$nnd
    if(want.which) nnwh[oX] <- if(is.sorted.Y) nnwcode else oY[nnwcode]
  }

  if(want.both) return(data.frame(dist=nndv, which=nnwh))
  return(if(want.dist) nndv else nnwh)
}

