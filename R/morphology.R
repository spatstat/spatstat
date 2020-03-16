#
#  morphology.R
#
#  dilation, erosion, opening, closing
#
#  generic functions
#  and methods for owin, psp, ppp
#
#  $Revision: 1.33 $   $Date: 2020/03/16 10:28:51 $
#

# ............ generic  ............................

erosion  <- function(w, r, ...) { UseMethod("erosion") }

dilation <- function(w, r, ...) { UseMethod("dilation") }

closing  <- function(w, r, ...) { UseMethod("closing") }

opening  <- function(w, r, ...) { UseMethod("opening") }

# ............ methods for class 'owin' ............................


# DELETED
# erode.owin <- function(...) {
#   .Deprecated("erosion.owin", package="spatstat")
#   erosion.owin(...)
# }

erosion.owin <- 
  function(w, r, shrink.frame=TRUE, ..., strict=FALSE, polygonal=NULL) {
  verifyclass(w, "owin")
  validradius(r, "erosion")
  if(r == 0 && !strict)
    return(w)

  xr <- w$xrange
  yr <- w$yrange
  
  if(2 * r >= max(diff(xr), diff(yr)))
    stop("erosion distance r too large for frame of window")

  # compute the dimensions of the eroded frame
  exr <- xr + c(r, -r)
  eyr <- yr + c(r, -r)
  ebox <- list(x=exr[c(1,2,2,1)], y=eyr[c(1,1,2,2)])

  ismask <- is.mask(w)
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
      if(is.mask(w))
        polygonal <- FALSE
    }
  }
  
  if(is.rectangle(w) && polygonal) {
    # result is a smaller rectangle
    if(shrink.frame) {
      return(owin(exr, eyr))  # type 'rectangle' 
    } else {
      return(owin(xr, yr, poly=ebox, check=FALSE)) # type 'polygonal'
    }
  }

  if(polygonal) {
    # compute polygonal region using polyclip package
    pnew <- polyclip::polyoffset(w$bdry, -r, jointype="round")
    # ensure correct polarity
    totarea <- sum(unlist(lapply(pnew, Area.xypolygon)))
    if(totarea < 0)
      pnew <- lapply(pnew, reverse.xypolygon)
    if(shrink.frame) {
      return(owin(poly=pnew, check=FALSE))
    } else {
      return(owin( xr,  yr, poly=pnew, check=FALSE))
    }
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

# DELETED
# dilate.owin <- function(...) {
#   .Deprecated("dilation.owin", package="spatstat")
#   dilation.owin(...)
# }

dilation.owin <- 
  function(w, r, ..., polygonal=NULL, tight=TRUE) {
  verifyclass(w, "owin")
  validradius(r, "dilation")
  
  if(r == 0)
    return(w)

  ismask <- is.mask(w)
  if(is.empty(w))
    return(w)

  # determine type of computation
  if(is.null(polygonal)) {
    polygonal <- !ismask
  } else stopifnot(is.logical(polygonal))
  
  if(polygonal) {
    # convert to polygonal 
    w <- as.polygonal(w)
    if(!is.polygonal(w))
      polygonal <- FALSE
  }
  
  # bounding frame
  bb <- if(tight) boundingbox(w) else as.rectangle(w)
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
    # compute polygonal region using polyclip package
    pnew <- polyclip::polyoffset(w$bdry, r, jointype="round")
    # ensure correct polarity
    totarea <- sum(unlist(lapply(pnew, Area.xypolygon)))
    if(totarea < 0)
      pnew <- lapply(pnew, reverse.xypolygon)
    # determine bounding frame, convert to owin
    if(tight) {
      out <- owin(poly=pnew, check=FALSE)
    } else {
      out <- owin(newbox$xrange, newbox$yrange, poly=pnew, check=FALSE)
    }
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
  wclose <- erosion.owin(wplus, r, strict=TRUE)
  b <- as.rectangle(w)
  wclose <- rebound.owin(wclose[b], b)
  return(wclose)
}

opening.owin <- function(w, r, ..., polygonal=NULL) {
  if(missing(r))
    stop("r is required")
  validradius(r, "opening")
  wminus <- erosion.owin(w, r, ..., polygonal=polygonal, shrink.frame=FALSE)
  if(is.empty(wminus))
    return(wminus)
  wopen <- dilation.owin(wminus, r, tight=FALSE)
  b <- as.rectangle(w)
  wopen <- rebound.owin(wopen[b], b)
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
  bb <- if(tight) boundingbox(x) else as.rectangle(x)
  newbox <- grow.rectangle(bb, r)
  
  # compute dilation
  if(!polygonal) {
    x <- rebound.psp(x, newbox)
    distant <- distmap(x, ...)
    dil <- levelset(distant, r, "<=")
    return(dil)
  } else if(spatstat.options("old.morpho.psp")) {
    # old code for polygonal case
    ends   <- x$ends
    angles <- angles.psp(x, directed=TRUE)
#    lengths <- lengths_psp(x)
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
      sausage <- owin(newbox$xrange, newbox$yrange, poly=list(x=xx, y=yy), check=FALSE)
      # add to set
      out <- union.owin(out, sausage, ...)
    }
    return(out)
  } else {
    # new code using 'polyclip' package
    # convert to list of list(x,y)
    ends   <- as.matrix(x$ends)
    n <- nrow(ends)
    plines <- vector(mode="list", length=n)
    for(i in 1:n) plines[[i]] <- list(x=ends[i, c("x0","x1")],
                                      y=ends[i, c("y0","y1")])
    # call
    pnew <- polyclip::polylineoffset(plines, r,
                                     jointype="round", endtype="openround")
    # ensure correct polarity
    totarea <- sum(unlist(lapply(pnew, Area.xypolygon)))
    if(totarea < 0)
      pnew <- lapply(pnew, reverse.xypolygon)
    # convert to owin object
    out <- if(tight) owin(poly=pnew, check=FALSE) else
            owin(newbox$xrange, newbox$yrange, poly=pnew, check=FALSE)
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
  bb <- if(tight) boundingbox(x) else as.rectangle(x)
  releps <- 1e-6
  newbox <- grow.rectangle(bb, r * (1+releps))

  # compute dilation
  if(!polygonal) {
    # compute pixel approximation
    Window(x) <- newbox
    distant <- distmap(x, ...)
    dil <- levelset(distant, r, "<=")
    return(dil)
  } else {
    # compute polygonal approximation
    # generate discs
    coo <- coords(x)
    nn <- npoints(x)
    balls <- vector(mode="list", length=nn)
    ball0 <- disc(r, c(0,0), ...)
    for(i in seq_len(nn))
      balls[[i]] <- shift(ball0, vec=coo[i,])
    class(balls) <- c("solist", class(balls))
    out <- union.owin(balls)
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

validradius <- local({

  validradius <- function(r, caller="morphological operator") {
  #  rname <- short.deparse(substitute(r))
    if(!is.numeric(r) || length(r) != 1)
      groan("radius r must be a single number", caller)
    if(r < 0)
      groan("radius r must be nonnegative", caller)
    return(TRUE)
  }

  groan <- function(whinge, caller) {
    stop(paste("for", paste(caller, ",", sep=""), whinge), call.=FALSE)
  }

  validradius
})
  
idorempty <- function(w, r, caller="morphological operator") {
  validradius(r, caller)
  if(r == 0)
    return(w)
  else
    return(emptywindow(w))
}
           
