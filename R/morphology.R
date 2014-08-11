#
#  morphology.R
#
#  dilation, erosion, opening, closing
#
#  generic functions
#  and methods for owin, psp, ppp
#
#  $Revision: 1.19 $   $Date: 2013/02/14 01:41:13 $
#

# ............ generic  ............................

erosion  <- function(w, r, ...) { UseMethod("erosion") }

dilation <- function(w, r, ...) { UseMethod("dilation") }

closing  <- function(w, r, ...) { UseMethod("closing") }

opening  <- function(w, r, ...) { UseMethod("opening") }

# ............ methods for class 'owin' ............................


erode.owin <- function(...) {
  .Deprecated("erosion.owin", package="spatstat")
  erosion.owin(...)
}

erosion.owin <- 
  function(w, r, shrink.frame=TRUE, ..., strict=FALSE, polygonal=NULL) {
  verifyclass(w, "owin")
  validradius(r, "erosion")
  if(r == 0 && !strict)
    return(w)
  
  if(2 * r >= max(diff(w$xrange), diff(w$yrange)))
    stop("erosion distance r too large for frame of window")

  # compute the dimensions of the eroded frame
  exr <- w$xrange + c(r, -r)
  eyr <- w$yrange + c(r, -r)
  ebox <- list(x=exr[c(1,2,2,1)], y=eyr[c(1,1,2,2)])

  ismask <- (w$type == "mask")
  if(is.empty(w))
    return(emptywindow(ebox))

  # determine type of computation
  if(is.null(polygonal))
    polygonal <- !ismask
  else {
    stopifnot(is.logical(polygonal))
    if(polygonal && ismask) {
      # try to convert
      w <- as.polygonal(w)
      if(w$type == "mask")
        polygonal <- FALSE
    }
  }
  
  if(w$type == "rectangle" && polygonal) {
    # result is a smaller rectangle
    if(shrink.frame)
      return(owin(exr, eyr))  # type 'rectangle' 
    else
      return(owin(w$xrange, w$yrange, poly=ebox)) # type 'polygonal'
  }

  if(polygonal) {
    # compute polygonal region
    cw <- complement.owin(w)
    dcw <- dilation.owin(cw, r, polygonal=TRUE, ...)
    cdcw <- complement.owin(dcw)
    wnew <- intersect.owin(cdcw, w)
    return(wnew)
  }
  
  # otherwise erode the window in pixel image form
  if(w$type == "mask") 
    wnew <- erodemask(w, r, strict=strict)
  else {
    D <- distmap(w, invert=TRUE, ...)
    wnew <- levelset(D, r, if(strict) ">" else ">=")
  }
        
  if(shrink.frame) {
    # trim off some rows & columns of pixel raster
    keepcol <- (wnew$xcol >= exr[1] & wnew$xcol <= exr[2])
    keeprow <- (wnew$yrow >= eyr[1] & wnew$yrow <= eyr[2])
    wnew$xcol <- wnew$xcol[keepcol]
    wnew$yrow <- wnew$yrow[keeprow]
    wnew$dim <- c(sum(keeprow), sum(keepcol))
    wnew$m <- wnew$m[keeprow, keepcol]
    wnew$xrange <- exr
    wnew$yrange <- eyr
  }

  return(wnew)
}	

dilate.owin <- function(...) {
  .Deprecated("dilation.owin", package="spatstat")
  dilation.owin(...)
}

dilation.owin <- 
  function(w, r, ..., polygonal=NULL, tight=TRUE) {
  verifyclass(w, "owin")
  validradius(r, "dilation")
  
  if(r == 0)
    return(w)

  ismask <- (w$type == "mask")
  if(is.empty(w))
    return(w)

  # determine type of computation
  if(is.null(polygonal))
    polygonal <- !ismask
  else {
    stopifnot(is.logical(polygonal))
    if(polygonal && ismask) {
      # try to convert
      w <- as.polygonal(w)
      if(w$type == "mask")
        polygonal <- FALSE
    }
  }
  
  
  # bounding frame
  bb <- if(tight) bounding.box(w) else as.rectangle(w)
  newbox <- grow.rectangle(bb, r)

  # compute dilation
  if(!polygonal) {
    # compute pixel approximation
    epsilon <- sqrt(w$xstep^2 + w$ystep^2)
    r <- max(r, epsilon)
    w <- rebound.owin(w, newbox)
    distant <- distmap(w, ...)
    dil <- levelset(distant, r, "<=")
    return(dil)
  } else {
    # compute polygonal approximation
    # extract individual edges
    w <- as.polygonal(w)
    edges <- as.psp(w)
    # dilate edges
    edgesplus <- dilation(edges, r, ..., polygonal=polygonal, tight=tight)
    # add interior of w
    out <- union.owin(w, edgesplus)
    return(out)
  }
}

closing.owin <- function(w, r, ..., polygonal=NULL) {
  if(missing(r))
    stop("r is required")
  validradius(r, "closing")
  wplus <- dilation.owin(w, r, ..., polygonal=polygonal, tight=FALSE)
  if(is.empty(wplus))
    return(wplus)
  if(wplus$type == "polygonal")
    wplus <- simplify.owin(wplus, r/10)
  wclose <- erosion.owin(wplus, r, strict=TRUE)
  wclose <- rebound.owin(wclose, as.rectangle(w))
  return(wclose)
}

opening.owin <- function(w, r, ..., polygonal=NULL) {
  if(missing(r))
    stop("r is required")
  validradius(r, "opening")
  wminus <- erosion.owin(w, r, ..., polygonal=polygonal, shrink.frame=FALSE)
  if(is.empty(wminus))
    return(wminus)
  if(wminus$type == "polygonal")
    wminus <- simplify.owin(wminus, r/10)
  wopen <- dilation.owin(wminus, r, tight=FALSE)
  wopen <- rebound.owin(wopen, as.rectangle(w))
  return(wopen)
}


border <- function(w, r, outside=FALSE, ...) {
  w <- as.owin(w)
  if(!outside) {
    e <- erosion(w, r, ...)
    b <- setminus.owin(w, e)
  } else {
    d <- dilation(w, r, ...)
    b <- setminus.owin(d, w)
  }
  return(b)
}

# ............ methods for class 'psp' ............................


dilation.psp <- function(w, r, ..., polygonal=TRUE, tight=TRUE) {
  verifyclass(w, "psp")
  x <- w
  validradius(r, "dilation")
  if(r == 0)
    return(w)

  if(is.empty(x))
    return(emptywindow(as.owin(w)))
  
  # bounding frame
  bb <- if(tight) bounding.box(x) else as.rectangle(x)
  newbox <- grow.rectangle(bb, r)
  
  # compute dilation
  if(!polygonal) {
    x <- rebound.psp(x, newbox)
    distant <- distmap(x, ...)
    dil <- levelset(distant, r, "<=")
    return(dil)
  } else {
    ends   <- x$ends
    angles <- angles.psp(x, directed=TRUE)
    lengths <- lengths.psp(x)
    out <- NULL
    # dilate individual segments
    halfcircle <- seq(from=0, to=pi, length.out=128)[-c(1,128)]
    for(i in seq_len(x$n)) {
      seg <- ends[i,]
      co <- cos(angles[i])
      si <- sin(angles[i])
      # draw sausage around i-th segment
      xx <- c(seg$x0, seg$x1) + r * si
      yy <- c(seg$y0, seg$y1) - r * co
      rightcircle <- angles[i] - pi/2 + halfcircle
      xx <- c(xx, seg$x1 + r * cos(rightcircle))
      yy <- c(yy, seg$y1 + r * sin(rightcircle))
      xx <- c(xx, c(seg$x1, seg$x0) - r * si)
      yy <- c(yy, c(seg$y1, seg$y0) + r * co)
      leftcircle <- angles[i] + pi/2 + halfcircle
      xx <- c(xx, seg$x0 + r * cos(leftcircle))
      yy <- c(yy, seg$y0 + r * sin(leftcircle))
      sausage <- owin(newbox$xrange, newbox$yrange, poly=list(x=xx, y=yy))
      # add to set
      out <- union.owin(out, sausage, ...)
    }
    return(out)
  }
}

closing.psp <- function(w, r, ..., polygonal=TRUE) {
  if(missing(r))
    stop("r is required")
  validradius(r, "closing")
  wplus <- dilation.psp(w, r, ..., polygonal=polygonal, tight=FALSE)
  if(is.empty(wplus))
    return(emptywindow(as.owin(w)))
  if(wplus$type == "polygonal")
    wplus <- simplify.owin(wplus, r/10)
  wclose <- erosion.owin(wplus, r, strict=TRUE)
  wclose <- rebound.owin(wclose, as.rectangle(w))
  return(wclose)
}

erosion.psp <- function(w, r, ...) {
  idorempty(w, r, "erosion")
}

opening.psp <- function(w, r, ...) {
  idorempty(w, r,"opening")
}

  
# ............ methods for class 'ppp' ............................

dilation.ppp <- function(w, r, ..., polygonal=TRUE, tight=TRUE) {
  verifyclass(w, "ppp")
  validradius(r, "dilation")
  x <- w
  
  if(r == 0)
    return(x)

  if(is.empty(w))
    return(emptywindow(as.owin(w)))
  
  # bounding frame
  bb <- if(tight) bounding.box(x) else as.rectangle(x)
    newbox <- grow.rectangle(bb, r)

  # compute dilation
  if(!polygonal) {
    # compute pixel approximation
    x <- rebound.ppp(x, newbox)
    distant <- distmap(x, ...)
    dil <- levelset(distant, r, "<=")
    return(dil)
  } else {
    # compute polygonal approximation
    # generate discs
    out <- NULL
    for(i in seq_len(x$n)) {
      balli <- disc(r, c(x$x[i], x$y[i]))
      out <- union.owin(out, balli)
    }
    return(out)
  }
}

closing.ppp <- function(w, r, ..., polygonal=TRUE) {
  if(missing(r))
    stop("r is required")
  validradius(r, "closing")
  if(is.empty(w) || w$n <= 3)
    return(emptywindow(as.owin(w)))
  # remove `isolated' points
  ok <- (nndist(w) <= 2 * r)
  if(sum(ok) <= 3)
    return(emptywindow(as.owin(w)))
  w <- w[ok]
  # dilate
  wplus <- dilation.ppp(w, r, ..., polygonal=polygonal, tight=FALSE)
  if(wplus$type == "polygonal")
    wplus <- simplify.owin(wplus, r/10)
  wclose <- erosion.owin(wplus, r, strict=TRUE)
  wclose <- rebound.owin(wclose, as.rectangle(w))
  return(wclose)
}

erosion.ppp <- function(w, r, ...) {
  idorempty(w, r, "erosion")
}

opening.ppp <- function(w, r, ...) {
  idorempty(w, r,"opening")
}

# ............ utilities ............................

validradius <- function(r, caller="morphological operator") {
  rname <- short.deparse(substitute(r))
  groan <- function(whinge, caller) {
    stop(paste("for", paste(caller, ",", sep=""), whinge), call.=FALSE)
  }
  if(!is.numeric(r) || length(r) != 1)
    groan("radius r must be a single number", caller)
  if(r < 0)
    groan("radius r must be nonnegative", caller)
  return(TRUE)
}
  
idorempty <- function(w, r, caller="morphological operator") {
  validradius(r, caller)
  if(r == 0)
    return(w)
  else
    return(emptywindow(w))
}
           
