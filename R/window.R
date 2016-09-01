#
#	window.S
#
#	A class 'owin' to define the "observation window"
#
#	$Revision: 4.175 $	$Date: 2016/09/01 05:52:38 $
#
#
#	A window may be either
#
#		- rectangular:
#                       a rectangle in R^2
#                       (with sides parallel to the coordinate axes)
#
#		- polygonal:
#			delineated by 0, 1 or more non-self-intersecting
#                       polygons, possibly including polygonal holes.
#	
#		- digital mask:
#			defined by a binary image
#			whose pixel values are TRUE wherever the pixel
#                       is inside the window
#
#	Any window is an object of class 'owin', 
#       containing at least the following entries:	
#
#		$type:	a string ("rectangle", "polygonal" or "mask")
#
#		$xrange   
#		$yrange
#			vectors of length 2 giving the real dimensions 
#			of the enclosing box
#               $units
#                       name of the unit of length
#
#	The 'rectangle' type has only these entries.
#
#       The 'polygonal' type has an additional entry
#
#               $bdry
#                       a list of polygons.
#                       Each entry bdry[[i]] determines a closed polygon.
#
#                       bdry[[i]] has components $x and $y which are
#                       the cartesian coordinates of the vertices of
#                       the i-th boundary polygon (without repetition of
#                       the first vertex, i.e. same convention as in the
#                       plotting function polygon().)
#
#
#	The 'mask' type has entries
#
#		$m		logical matrix
#		$dim		its dimension array
#		$xstep,ystep	x and y dimensions of a pixel
#		$xcol	        vector of x values for each column
#               $yrow           vector of y values for each row
#	
#	(the row index corresponds to increasing y coordinate; 
#	 the column index "   "     "   "  "  "  x "   "    ".)
#
#
#-----------------------------------------------------------------------------
#

.Spatstat.Image.Warning <-
  c("Row index corresponds to increasing y coordinate; column to increasing x",
    "Transpose matrices to get the standard presentation in R",
    "Example: image(result$xcol,result$yrow,t(result$d))")

owin <- local({

  isxy <- function(x) { (is.matrix(x) || is.data.frame(x)) && ncol(x) == 2 }
  asxy <- function(xy) { list(x=xy[,1], y=xy[,2]) }

  owin <- function(xrange=c(0,1), yrange=c(0,1),
                 ..., poly=NULL, mask=NULL, unitname=NULL, xy=NULL) {

  # trap a common abuse of syntax
  if(nargs() == 1 && !missing(xrange) && is.owin(xrange))
    return(xrange)
  
  unitname <- as.units(unitname)

  ## Exterminate ambiguities
  poly.given <- !is.null(poly)
  mask.given <- !is.null(mask)
  if(poly.given  && mask.given)
     stop("Ambiguous -- both polygonal boundary and digital mask supplied")

  if(!is.null(xy) && !mask.given)
    warning("Argument xy ignored: it is only applicable when a mask is given")
     
  if(missing(xrange) != missing(yrange))
    stop("If one of xrange, yrange is specified then both must be.")

  # convert data frames to vanilla lists
  if(poly.given) {
    if(is.data.frame(poly))
      poly <- as.list(poly)
    else if(is.list(poly) && any(unlist(lapply(poly, is.data.frame))))
      poly <- lapply(poly, as.list)
  }
  
  ## Hidden options controlling how much checking is performed
  check <- resolve.1.default(list(check=TRUE), list(...))
  calculate <- resolve.1.default(list(calculate=check), list(...))
  strict <- resolve.1.default(list(strict=spatstat.options("checkpolygons")),
                              list(...))
  fix <- resolve.1.default(list(fix=spatstat.options("fixpolygons")),
                           list(...))

  if(!poly.given && !mask.given) {
    ######### rectangle #################
    if(check) {
      if(!is.vector(xrange) || length(xrange) != 2 || xrange[2] < xrange[1])
        stop("xrange should be a vector of length 2 giving (xmin, xmax)")
      if(!is.vector(yrange) || length(yrange) != 2 || yrange[2] < yrange[1])
        stop("yrange should be a vector of length 2 giving (ymin, ymax)")
    }
    w <- list(type="rectangle", xrange=xrange, yrange=yrange, units=unitname)
    class(w) <- "owin"
    return(w)
  } else if(poly.given) {
    ######### polygonal boundary ########
    #
    if(length(poly) == 0) {
      # empty polygon
      if(check) {
        if(!is.vector(xrange) || length(xrange) != 2 || xrange[2] < xrange[1])
          stop("xrange should be a vector of length 2 giving (xmin, xmax)")
        if(!is.vector(yrange) || length(yrange) != 2 || yrange[2] < yrange[1])
          stop("yrange should be a vector of length 2 giving (ymin, ymax)")
      }
      w <- list(type="polygonal", xrange=xrange, yrange=yrange,
                bdry=list(), units=unitname)
      class(w) <- "owin"
      return(w)
    }
    # convert matrix or data frame to list(x,y)
    if(isxy(poly)) {
      poly <- asxy(poly)
    } else if(is.list(poly) && all(unlist(lapply(poly, isxy)))) {
      poly <- lapply(poly, asxy)
    }
    # nonempty polygon  
    # test whether it's a single polygon or multiple polygons
    if(verify.xypolygon(poly, fatal=FALSE))
      psingle <- TRUE
    else if(all(unlist(lapply(poly, verify.xypolygon, fatal=FALSE))))
      psingle <- FALSE
    else
      stop("poly must be either a list(x,y) or a list of list(x,y)")

    w.area <- NULL
    
    if(psingle) {
      # single boundary polygon
      bdry <- list(poly)
      if(check || calculate) {
        w.area <- Area.xypolygon(poly)
        if(w.area < 0)
          stop(paste("Area of polygon is negative -",
                     "maybe traversed in wrong direction?"))
      }
    } else {
      # multiple boundary polygons
      bdry <- poly
      if(check || calculate) {
        w.area <- unlist(lapply(poly, Area.xypolygon))
        if(sum(w.area) < 0)
          stop(paste("Area of window is negative;\n",
                     "check that all polygons were traversed",
                     "in the right direction"))
      }
    }

    actual.xrange <- range(unlist(lapply(bdry, getElement, name="x")))
    if(missing(xrange))
      xrange <- actual.xrange
    else if(check) {
      if(!is.vector(xrange) || length(xrange) != 2 || xrange[2] <= xrange[1])
        stop("xrange should be a vector of length 2 giving (xmin, xmax)")
      if(!all(xrange == range(c(xrange, actual.xrange))))
        stop("polygon's x coordinates outside xrange")
    }
    
    actual.yrange <- range(unlist(lapply(bdry, getElement, name="y")))
    if(missing(yrange))
      yrange <- actual.yrange
    else if(check) {
      if(!is.vector(yrange) || length(yrange) != 2 || yrange[2] <= yrange[1])
        stop("yrange should be a vector of length 2 giving (ymin, ymax)")
      if(!all(yrange == range(c(yrange, actual.yrange))))
      stop("polygon's y coordinates outside yrange")
    }

    if(!is.null(w.area)) {
      # tack on area and hole data
      holes <- (w.area < 0)
      for(i in seq_along(bdry)) 
        bdry[[i]] <- append(bdry[[i]], list(area=w.area[i], hole=holes[i]))
    }
    
    w <- list(type="polygonal",
              xrange=xrange, yrange=yrange, bdry=bdry, units=unitname)
    class(w) <- "owin"

    if(check && strict) { 
      ## strict checks on geometry (self-intersection etc)
      ok <- owinpolycheck(w)
      if(!ok) {
        errors <- attr(ok, "err")
        stop(paste("Polygon data contain", commasep(errors)))
      }
    }
    if(check && fix) {
      if(length(bdry) == 1 &&
         length(bx <- bdry[[1]]$x) == 4 &&
         length(unique(bx)) == 2 &&
         length(unique(bdry[[1]]$y)) == 2) {
        ## it's really a rectangle
        if(Area.xypolygon(bdry[[1]]) < 0)
          w$bdry <- lapply(bdry, reverse.xypolygon)
      } else {
        ## repair polygon data by invoking polyclip
        ##        to intersect polygon with larger-than-bounding rectangle
        ##        (Streamlined version of intersect.owin)
        ww <- lapply(bdry, reverse.xypolygon)
        xrplus <- mean(xrange) + c(-1,1) * diff(xrange)
        yrplus <- mean(yrange) + c(-1,1) * diff(yrange)
        bignum <- (.Machine$integer.max^2)/2
        epsclip <- max(diff(xrange), diff(yrange))/bignum
        rr <- list(list(x=xrplus[c(1,2,2,1)], y=yrplus[c(2,2,1,1)]))
        bb <- polyclip::polyclip(ww, rr, "intersection",
                                 fillA="nonzero", fillB="nonzero", eps=epsclip)
        ## ensure correct polarity
        totarea <- sum(unlist(lapply(bb, Area.xypolygon)))
        if(totarea < 0)
          bb <- lapply(bb, reverse.xypolygon)
        w$bdry <- bb
      }
    }
    return(w)
  } else if(mask.given) {
    ######### digital mask #####################

    if(is.data.frame(mask) &&
       ncol(mask) %in% c(2,3) &&
       sum(sapply(mask, is.numeric)) == 2) {
      # data frame with 2 columns of coordinates
      return(as.owin(W=mask, xy=xy))
    }
      
    if(!is.matrix(mask))
      stop(paste(sQuote("mask"), "must be a matrix"))
    if(!is.logical(mask))
      stop(paste("The entries of", sQuote("mask"), "must be logical"))
    
    nc <- ncol(mask)
    nr <- nrow(mask)

    if(!is.null(xy)) {
      # pixel coordinates given explicitly
      # validate dimensions
      if(!is.list(xy) || !checkfields(xy, c("x","y")))
        stop("xy should be a list with entries x and y")
      xcol <- xy$x
      yrow <- xy$y
      if(length(xcol) != nc)
        stop(paste("length of xy$x =", length(xcol),
                   "!=", nc, "= number of columns of mask"))
      if(length(yrow) != nr)
        stop(paste("length of xy$y =", length(yrow),
                   "!=", nr, "= number of rows of mask"))
      # x and y should be evenly spaced
      if(!evenly.spaced(xcol))
        stop("xy$x is not evenly spaced")
      if(!evenly.spaced(yrow))
        stop("xy$y is not evenly spaced")
      # determine other parameters
      xstep <- diff(xcol)[1]
      ystep <- diff(yrow)[1]
      if(missing(xrange) && missing(yrange)) {
        xrange <- range(xcol) + c(-1,1) * xstep/2
        yrange <- range(yrow) + c(-1,1) * ystep/2
      }
    } else {
      # determine pixel coordinates from xrange, yrange
      if(missing(xrange) && missing(yrange)) {
        # take pixels to be 1 x 1 unit
        xrange <- c(0,nc)
        yrange <- c(0,nr)
      } else if(check) {
        if(!is.vector(xrange) || length(xrange) != 2 || xrange[2] <= xrange[1])
          stop("xrange should be a vector of length 2 giving (xmin, xmax)")
        if(!is.vector(yrange) || length(yrange) != 2 || yrange[2] <= yrange[1])
          stop("yrange should be a vector of length 2 giving (ymin, ymax)")
      }
      xstep <- diff(xrange)/nc
      ystep <- diff(yrange)/nr
      xcol  <- seq(from=xrange[1]+xstep/2, to=xrange[2]-xstep/2, length.out=nc)
      yrow  <- seq(from=yrange[1]+ystep/2, to=yrange[2]-ystep/2, length.out=nr)
    }

    out <- list(type     = "mask",
                xrange   = xrange,
                yrange   = yrange,
                dim      = c(nr, nc),
                xstep    = xstep,
                ystep    = ystep,
                warnings = .Spatstat.Image.Warning,
                xcol    = xcol, 
                yrow    = yrow,
                m       = mask,
                units   = unitname)
    class(out) <- "owin"
    return(out)
  }
  # never reached
  NULL
}

  owin
})

#
#-----------------------------------------------------------------------------
#

is.owin <- function(x) { inherits(x, "owin") }

#
#-----------------------------------------------------------------------------
#

as.owin <- function(W, ..., fatal=TRUE) {
  UseMethod("as.owin")
}

as.owin.owin <- function(W, ..., fatal=TRUE) {
  if(verifyclass(W, "owin", fatal=fatal)) 
    return(owin(W$xrange, W$yrange, poly=W$bdry, mask=W$m, unitname=unitname(W), check=FALSE))
  else
    return(NULL)
}

as.owin.ppp <- function(W, ..., fatal=TRUE) {
  if(verifyclass(W, "ppp", fatal=fatal))
    return(W$window)
  else
    return(NULL)
}

as.owin.quad <- function(W, ..., fatal=TRUE) {
  if(verifyclass(W, "quad", fatal=fatal))
    return(W$data$window)
  else
    return(NULL)
}

as.owin.im <- function(W, ..., fatal=TRUE) {
  if(!verifyclass(W, "im", fatal=fatal))
    return(NULL)
  out <- list(type     = "mask",
              xrange   = W$xrange,
              yrange   = W$yrange,
              dim      = W$dim,
              xstep    = W$xstep,
              ystep    = W$ystep,
              warnings = .Spatstat.Image.Warning,
              xcol    = W$xcol,
              yrow    = W$yrow,
              m       = !is.na(W$v),
              units   = unitname(W))
  class(out) <- "owin"
  return(out)
}

as.owin.psp <- function(W, ..., fatal=TRUE) {
  if(!verifyclass(W, "psp", fatal=fatal))
    return(NULL)
  return(W$window)
}

as.owin.tess <- function(W, ..., fatal=TRUE) {
  if(!verifyclass(W, "tess", fatal=fatal))
    return(NULL)
  return(W$window)
}

as.owin.data.frame <- function(W, ..., step, fatal=TRUE) {
  if(!verifyclass(W, "data.frame", fatal=fatal))
    return(NULL)
  if(missing(step)) {
    xstep <- ystep <- NULL
  } else {
    step <- ensure2vector(step)
    xstep <- step[1]
    ystep <- step[2]
  }
  if(!(ncol(W) %in% c(2,3))) {
    whinge <- "need exactly 2 or 3 columns of data"
    if(fatal) stop(whinge)
    warning(whinge)
    return(NULL)
  }
  if(twocol <- (ncol(W) == 2)) {
    # assume data is a list of TRUE pixels
    W <- cbind(W, TRUE)
  } 
  mch <- matchNameOrPosition(c("x", "y", "z"), names(W))
  ix <- mch[1]
  iy <- mch[2]
  iz <- mch[3]
  df <- data.frame(x=W[,ix], y=W[,iy], z=as.logical(W[,iz]))
  with(df, {
    xx <- sort(unique(x))
    yy <- sort(unique(y))
    jj <- match(x, xx)
    ii <- match(y, yy)
    ## make logical matrix (for incomplete x, y sequence)
    ok <- checkbigmatrix(length(xx), length(yy), fatal=fatal)
    if(!ok) return(NULL)
    mm <- matrix(FALSE, length(yy), length(xx))
    mm[cbind(ii,jj)] <- z
    ## ensure xx and yy are complete equally-spaced sequences
    fx <- fillseq(xx, step=xstep)
    fy <- fillseq(yy, step=ystep)
    xcol <- fx[[1]]
    yrow <- fy[[1]]
    ## trap very large matrices
    ok <- checkbigmatrix(length(xcol), length(yrow), fatal=fatal)
    if(!ok) return(NULL)
    ## mapping from xx to xcol, yy to yrow
    jjj <- fx[[2]]
    iii <- fy[[2]]
    ## make logical matrix for full sequence
    m <- matrix(FALSE, length(yrow), length(xcol))
    m[iii,jjj] <- mm
    ## make binary mask
    out <- owin(mask=m, xy=list(x=xcol, y=yrow))
    ## warn if area fraction is small: may be a misuse of as.owin
    if(twocol) {
      pcarea <- 100 * nrow(df)/prod(dim(m))
      if(pcarea < 1) 
       warning(paste("Window occupies only",
                     paste0(signif(pcarea, 2), "%"),
                     "of frame area. Did you mean owin(poly=df) ?"),
               call.=FALSE)
    }
    return(out)
  })
}

as.owin.default <- function(W, ..., fatal=TRUE) {
  ## Tries to interpret data as an object of class 'owin'
  ## W may be
  ##	a structure with entries xrange, yrange
  ##	a four-element vector (interpreted xmin, xmax, ymin, ymax)
  ##	a structure with entries xl, xu, yl, yu
  ##	an object with attribute "bbox"

  if(checkfields(W, c("xrange", "yrange"))) {
    Z <- owin(W$xrange, W$yrange)
    return(Z)
  } else if(is.vector(W) && is.numeric(W) && length(W) == 4) {
    Z <- owin(W[1:2], W[3:4])
    return(Z)
  } else if(checkfields(W, c("xl", "xu", "yl", "yu"))) {
    W <- as.list(W)
    Z <- owin(c(W$xl, W$xu),c(W$yl, W$yu))
    return(Z)
  } else if(checkfields(W, c("x", "y", "area"))
            && checkfields(W$area, c("xl", "xu", "yl", "yu"))) {
    V <- as.list(W$area)
    Z <- owin(c(V$xl, V$xu),c(V$yl, V$yu))
    return(Z)
  } else if(!is.null(Z <- attr(W, "bbox"))) {
    return(as.owin(Z, ..., fatal=fatal))
  } else if(fatal)
    stop("Can't interpret W as a window")
  else return(NULL)
}		

#
#-----------------------------------------------------------------------------
#
#
Frame <- function(X) { UseMethod("Frame") }

"Frame<-" <- function(X, value) { UseMethod("Frame<-") }

Frame.default <- function(X) { as.rectangle(X) }

## .........................................................

as.rectangle <- function(w, ...) {
  if(inherits(w, "owin"))
    return(owin(w$xrange, w$yrange, unitname=unitname(w)))
  else if(inherits(w, "im"))
    return(owin(w$xrange, w$yrange, unitname=unitname(w)))
  else if(inherits(w, "layered")) 
    return(do.call(boundingbox, unname(lapply(w, as.rectangle))))
  else {
    w <- as.owin(w, ...)
    return(owin(w$xrange, w$yrange, unitname=unitname(w)))
  }
}

#
#-----------------------------------------------------------------------------
#
as.mask <- function(w, eps=NULL, dimyx=NULL, xy=NULL) {
#	eps:		   grid mesh (pixel) size
#	dimyx:		   dimensions of pixel raster
#       xy:                coordinates of pixel raster
  nonamedargs <- is.null(eps) && is.null(dimyx) && is.null(xy)
  uname <- as.units(NULL)
  if(!missing(w) && !is.null(w)) {
    if(is.data.frame(w)) return(owin(mask=w, xy=xy))
    if(is.matrix(w)) {
      w <- as.data.frame(w)
      colnames(w) <- c("x", "y")
      return(owin(mask=w, xy=xy))
    }
    w <- as.owin(w)
    uname <- unitname(w)
  } else {
    if(is.null(xy)) 
      stop("If w is missing, xy is required")
  }
  # If it's already a mask, and no other arguments specified,
  # just return it.
  if(!missing(w) && w$type == "mask" && nonamedargs)
    return(w)
  
##########################
#  First determine pixel coordinates
##########################
  if(is.null(xy)) {
# Pixel coordinates to be computed from other dimensions
# First determine row & column dimensions
    if(!is.null(dimyx)) {
      dimyx <- ensure2vector(dimyx)
      nr <- dimyx[1]
      nc <- dimyx[2]
    } else {
    # use pixel size 'eps'
      if(!is.null(eps)) {
        eps <- ensure2vector(eps)
        nc <- diff(w$xrange)/eps[1]
        nr <- diff(w$yrange)/eps[2]
        if(nr < 1 || nc < 1)
          warning("pixel size parameter eps > size of window")
        nr <- ceiling(nr)
        nc <- ceiling(nc)
      } else {
    # use spatstat defaults
        np <- spatstat.options("npixel")
        if(length(np) == 1)
          nr <- nc <- np[1]
        else {
          nr <- np[2]  
          nc <- np[1]
        }
      }
    }
    if((mpix <- (nr * nc)/1048576) >= 10) {
      whinge <- paste("Creating",
                      articlebeforenumber(mpix),
                      paste0(round(mpix, 1), "-megapixel"),
                      "window mask")
      message(whinge)
      warning(whinge, call.=FALSE)
    }
    # Initialise mask with all entries TRUE
    rasta <- owin(w$xrange, w$yrange, mask=matrix(TRUE, nr, nc))
  } else {
# 
# Pixel coordinates given explicitly:
#    xy is an image, a mask, or a list(x,y)
#
    if(is.im(xy)) {
      rasta <- as.owin(xy)
      rasta$m[] <- TRUE
    } else if(is.owin(xy)) {
      if(xy$type != "mask")
        stop("argument xy does not contain raster coordinates.")
      rasta <- xy
      rasta$m[] <- TRUE
    } else {
      if(!checkfields(xy, c("x", "y")))
        stop(paste(sQuote("xy"),
                   "should be a list containing two vectors x and y"))
      x <- sort(unique(xy$x))
      y <- sort(unique(xy$y))
      # derive other parameters
      nr <- length(y)
      nc <- length(x)
      # check size
      if((mpix <- (nr * nc)/1048576) >= 10) {
        whinge <- paste("Creating",
                        articlebeforenumber(mpix),
                        paste0(round(mpix, 1), "-megapixel"),
                        "window mask")
        message(whinge)
        warning(whinge, call.=FALSE)
      }
      # x and y pixel sizes
      dx <- diff(x)
      if(diff(range(dx)) > 0.01 * mean(dx))
        stop("x coordinates must be evenly spaced")
      xstep <- mean(dx)
      dy <- diff(y)
      if(diff(range(dy)) > 0.01 * mean(dy))
        stop("y coordinates must be evenly spaced")
      ystep <- mean(dy)
      xr <- range(x)
      yr <- range(y)
      xrange <-  xr + xstep * c(-1,1)/2
      yrange <-  yr + ystep * c(-1,1)/2
      # initialise mask with all entries TRUE
      rasta <- list(type     = "mask",
                    xrange   = xrange,
                    yrange   = yrange,
                    dim      = c(nr, nc),
                    xstep    = xstep,
                    ystep    = ystep,
                    warnings = .Spatstat.Image.Warning,
                    xcol    = seq(from=xr[1], to=xr[2], length.out=nc),
                    yrow    = seq(from=yr[1], to=yr[2], length.out=nr),
                    m       = matrix(TRUE, nr, nc),
                    units   = uname)
      class(rasta) <- "owin"
    }
    if(missing(w)) {
      # No more window information
      out <- rasta
      if(!(identical(x, xy$x) && identical(y, xy$y))) {
        ## xy is an enumeration of the TRUE pixels
        out$m[] <- FALSE
        ij <- cbind(i=match(xy$y, y), j=match(xy$x, x))
        out$m[ij] <- TRUE
      }
      return(out)
    }
  }

################################  
# Second, mask pixel raster with existing window
################################  
  switch(w$type,
         rectangle = {
           out <- rasta
           if(!all(w$xrange == rasta$xrange)
              || !all(w$yrange == rasta$yrange)) {
             xcol <- rasta$xcol
             yrow <- rasta$yrow
             wx <- w$xrange
             wy <- w$yrange
             badrow <- which(yrow > wy[2] | yrow < wy[1])
             badcol <- which(xcol > wx[2] | xcol < wx[1])
             out$m[badrow , ] <- FALSE
             out$m[ , badcol] <- FALSE
           }
         },
         mask = {
           # resample existing mask on new raster
           out <- rastersample(w, rasta)
         },
         polygonal = {
           # use C code
           out <- owinpoly2mask(w, rasta, FALSE)
         })

  unitname(out) <- uname
  return(out)
}

as.matrix.owin <- function(x, ...) {
  m <- as.mask(x, ...)
  return(m$m)
}

#
#
#-----------------------------------------------------------------------------
#
as.polygonal <- function(W, repair=FALSE) {
  verifyclass(W, "owin")
  switch(W$type,
         rectangle = {
           xr <- W$xrange
           yr <- W$yrange
           return(owin(xr, yr, poly=list(x=xr[c(1,2,2,1)],y=yr[c(1,1,2,2)]),
                       unitname=unitname(W),
                       check=FALSE))
         },
         polygonal = {
           if(repair)
             W <- owin(poly=W$bdry, unitname=unitname(W))
           return(W)
         },
         mask = {
           # This could take a while
           M <- W$m
           nr <- nrow(M)
           notM <- !M
           out <- NULL
           xcol <- W$xcol
           yrow <- W$yrow
           xbracket <- 1.1 * c(-1,1) * W$xstep/2
           ybracket <- 1.1 * c(-1,1) * W$ystep/2
           # identify runs of TRUE entries in each column
           start <- M & rbind(TRUE, notM[-nr, ])
           finish <- M & rbind(notM[-1, ], TRUE)
           for(j in 1:ncol(M)) {
             xj <- xcol[j]
             # identify start and end positions in column j
             starts <- which(start[,j])
             finishes <- which(finish[,j])
             ns <- length(starts)
             nf <- length(finishes)
             if(ns != nf)
               stop(paste("Internal error: length(starts)=", ns,
                          ", length(finishes)=", nf))
             if(ns > 0) {
               for(k in 1:ns) {
                 yfrom <- yrow[starts[k]]
                 yto   <- yrow[finishes[k]]
                 yk <- sort(c(yfrom,yto))
                 # make rectangle
                 recto <- owin(xj+xbracket,yk+ybracket)
                 # add to result
                 out <- union.owin(out, recto)
               }
               unitname(out) <- unitname(W)
             }
           }
           return(out)
         }
         )
}

#
# ----------------------------------------------------------------------

is.polygonal <- function(w) {
  return(inherits(w, "owin") && (w$type == "polygonal"))
}

is.rectangle <- function(w) {
  return(inherits(w, "owin") && (w$type == "rectangle"))
}

is.mask <- function(w) {
  return(inherits(w, "owin") && (w$type == "mask"))
}

validate.mask <- function(w, fatal=TRUE) {
  verifyclass(w, "owin", fatal=fatal)
  if(w$type == "mask")
    return(TRUE)
  if(fatal)
      stop(paste(short.deparse(substitute(w)), "is not a binary mask"))
  else {
      warning(paste(short.deparse(substitute(w)), "is not a binary mask"))
      return(FALSE)
  }
}

dim.owin <- function(x) { return(x$dim) } ## NULL unless it's a mask

## internal use only:

rasterx.mask <- function(w, drop=FALSE) {
  validate.mask(w)
  x <- w$xcol[col(w)]
  x <- if(drop) x[w$m, drop=TRUE] else array(x, dim=w$dim)
  return(x)
}

rastery.mask <- function(w, drop=FALSE) {
  validate.mask(w)
  y <- w$yrow[row(w)]
  y <- if(drop) y[w$m, drop=TRUE] else array(y, dim=w$dim)
  return(y)
}

rasterxy.mask <- function(w, drop=FALSE) {
  validate.mask(w)
  x <- w$xcol[col(w)]
  y <- w$yrow[row(w)]
  if(drop) {
    m <- w$m
    x <- x[m, drop=TRUE] 
    y <- y[m, drop=TRUE]
  }
  return(list(x=as.numeric(x),
              y=as.numeric(y)))
}


nearest.raster.point <- function(x,y,w, indices=TRUE) {
  stopifnot(is.mask(w) || is.im(w))
  nr <- w$dim[1]
  nc <- w$dim[2]
  if(length(x) == 0) {
    cc <- rr <- integer(0)
  } else {
    cc <- 1 + round((x - w$xcol[1])/w$xstep)
    rr <- 1 + round((y - w$yrow[1])/w$ystep)
    cc <- pmax.int(1,pmin.int(cc, nc))
    rr <- pmax.int(1,pmin.int(rr, nr))
  }
  if(indices) 
    return(list(row=rr, col=cc))
  else
    return(list(x=w$xcol[cc], y=w$yrow[rr]))
}

mask2df <- function(w) {
  stopifnot(is.owin(w) && w$type == "mask")
  xx <- raster.x(w)
  yy <- raster.y(w)
  ok <- w$m
  xx <- as.vector(xx[ok])
  yy <- as.vector(yy[ok])
  return(data.frame(x=xx, y=yy))
}

#------------------------------------------------------------------
		
complement.owin <- function(w, frame=as.rectangle(w)) {
  w <- as.owin(w)

  if(reframe <- !missing(frame)) {
    verifyclass(frame, "owin")
    w <- rebound.owin(w, frame)
    # if w was a rectangle, it's now polygonal
  }

  switch(w$type,
         mask = {
           w$m <- !(w$m)
         },
         rectangle = {
           # return empty window
           return(emptywindow(w))
         },
         polygonal = {
           bdry <- w$bdry
           if(length(bdry) == 0) {
             # w is empty
             return(frame)
           }
           # bounding box, in anticlockwise order
           box <- list(x=w$xrange[c(1,2,2,1)],
                       y=w$yrange[c(1,1,2,2)])
           boxarea <- Area.xypolygon(box)
                 
           # first check whether one of the current boundary polygons
           # is the bounding box itself (with + sign)
           if(reframe)
             is.box <- rep.int(FALSE, length(bdry))
           else {
             nvert <- lengths(lapply(bdry, getElement, name="x"))
             areas <- sapply(bdry, Area.xypolygon)
             boxarea.mineps <- boxarea * (0.99999)
             is.box <- (nvert == 4 & areas >= boxarea.mineps)
             if(sum(is.box) > 1)
               stop("Internal error: multiple copies of bounding box")
             if(all(is.box)) {
               return(emptywindow(box))
             }
           }
                 
           # if box is present (with + sign), remove it
           if(any(is.box))
             bdry <- bdry[!is.box]
                 
           # now reverse the direction of each polygon
           bdry <- lapply(bdry, reverse.xypolygon, adjust=TRUE)

           # if box was absent, add it
           if(!any(is.box))
             bdry <- c(bdry, list(box))   # sic

           # put back into w
           w$bdry <- bdry
         },
         stop("unrecognised window type", w$type)
         )
  return(w)
}

#-----------------------------------------------------------

inside.owin <- function(x, y, w) {
  # test whether (x,y) is inside window w
  # x, y may be vectors 

  if(missing(y) && all(c("x", "y") %in% names(x)))
    return(inside.owin(x$x, x$y, w))

  w <- as.owin(w)

  if(length(x)==0)
    return(logical(0))
  
  # test whether inside bounding rectangle
  xr <- w$xrange
  yr <- w$yrange
  eps <- sqrt(.Machine$double.eps)
  frameok <- (x >= xr[1] - eps) & (x <= xr[2] + eps) & 
             (y >= yr[1] - eps) & (y <= yr[2] + eps)
 
  if(!any(frameok))  # all points OUTSIDE window - no further work needed
    return(frameok)

  ok <- frameok
  switch(w$type,
         rectangle = {
           return(ok)
         },
         polygonal = {
           ## check scale
           framesize <- max(diff(xr), diff(yr))
           if(framesize > 1e6 || framesize < 1e-6) {
             ## rescale to avoid numerical overflow
             scalefac <- framesize/100
             w <- rescale(w, scalefac)
             x <- x/scalefac
             y <- y/scalefac
           }
           xy <- list(x=x,y=y)
           bdry <- w$bdry
           total <- numeric(length(x))
           on.bdry <- rep.int(FALSE, length(x))
           for(i in seq_along(bdry)) {
             score <- inside.xypolygon(xy, bdry[[i]], test01=FALSE)
             total <- total + score
             on.bdry <- on.bdry | attr(score, "on.boundary")
           }
           # any points identified as belonging to the boundary get score 1
           total[on.bdry] <- 1
           # check for sanity now..
           uhoh <- (total * (1-total) != 0)
           if(any(uhoh)) {
             nuh <- sum(uhoh)
             warning(paste("point-in-polygon test had difficulty with",
                           nuh,
                           ngettext(nuh, "point", "points"),
                           "(total score not 0 or 1)"),
                     call.=FALSE)
             total[uhoh] <- 0
           }
           return(ok & (total != 0))
         },
         mask = {
           # consider only those points which are inside the frame
           xf <- x[frameok]
           yf <- y[frameok]
           # map locations to raster (row,col) coordinates
           loc <- nearest.raster.point(xf,yf,w)
           # look up mask values
           okf <- (w$m)[cbind(loc$row, loc$col)]
           # insert into 'ok' vector
           ok[frameok] <- okf
           return(ok)
         },
         stop("unrecognised window type", w$type)
         )
}

#-------------------------------------------------------------------------
  
print.owin <- function(x, ..., prefix="window: ") {
  verifyclass(x, "owin")
  unitinfo <- summary(unitname(x))
  switch(x$type,
         rectangle={
           rectname <- paste0(prefix, "rectangle =")
         },
         polygonal={
           splat(paste0(prefix, "polygonal boundary"))
           if(length(x$bdry) == 0)
             splat("window is empty")
           rectname <- "enclosing rectangle:"
         },
         mask={
           splat(paste0(prefix, "binary image mask"))
           di <- x$dim
           splat(di[1], "x", di[2], "pixel array (ny, nx)")
           rectname <- "enclosing rectangle:"
         }
         )
  splat(rectname,
        prange(zapsmall(x$xrange)),
        "x",
        prange(zapsmall(x$yrange)),
        unitinfo$plural,
        unitinfo$explain)
  invisible(NULL)
}

summary.owin <- function(object, ...) {
  verifyclass(object, "owin")
  result <- list(xrange=object$xrange,
                 yrange=object$yrange,
                 type=object$type,
                 area=area(object),
                 units=unitname(object))
  result$areafraction <- with(result, area/(diff(xrange) * diff(yrange)))
  switch(object$type,
         rectangle={
         },
         polygonal={
           poly <- object$bdry
           result$npoly <- npoly <- length(poly)
           if(npoly == 0) {
             result$areas <- result$nvertices <- numeric(0)
           } else if(npoly == 1) {
             result$areas <- Area.xypolygon(poly[[1]])
             result$nvertices <- length(poly[[1]]$x)
           } else {
             result$areas <- unlist(lapply(poly, Area.xypolygon))
             result$nvertices <- lengths(lapply(poly, getElement, name="x"))
           }
           result$nhole <- sum(result$areas < 0)
         },
         mask={
           result$npixels <- object$dim
           result$xstep <- object$xstep
           result$ystep <- object$ystep
         }
         )
  class(result) <- "summary.owin"
  result
}

print.summary.owin <- function(x, ...) {
  verifyclass(x, "summary.owin")
  unitinfo <- summary(x$units)
  pluralunits <- unitinfo$plural
  singularunits <- unitinfo$singular
  switch(x$type,
         rectangle={
           rectname <- "Window: rectangle ="
         },
         polygonal={
           np <- x$npoly
           splat("Window: polygonal boundary")
           if(np == 0) {
             splat("window is empty")
           } else if(np == 1) {
             splat("single connected closed polygon with",
                   x$nvertices, "vertices")
           } else {
             nh <- x$nhole
             holy <- if(nh == 0) "(no holes)" else
                     if(nh == 1) "(1 hole)" else
                     paren(paste(nh, "holes"))
             splat(np, "separate polygons", holy)
             if(np > 0)
               print(data.frame(vertices=x$nvertices,
                                area=signif(x$areas, 6),
                                relative.area=signif(x$areas/x$area,3),
                                row.names=paste("polygon",
                                  1:np,
                                  ifelse(x$areas < 0, "(hole)", "")
                                  )))
           }
           rectname <- "enclosing rectangle:"
         },
         mask={
           splat("binary image mask")
           di <- x$npixels
           splat(di[1], "x", di[2], "pixel array (ny, nx)")
           splat("pixel size:",
                 signif(x$xstep,3), "by", signif(x$ystep,3),
                 pluralunits)
           rectname <- "enclosing rectangle:"
         }
         )
  splat(rectname,
        prange(zapsmall(x$xrange)),
        "x",
        prange(zapsmall(x$yrange)),
        pluralunits)
  Area <- signif(x$area, 6)
  splat("Window area =", Area, "square",
        if(Area == 1) singularunits else pluralunits)
  if(!is.null(ledge <- unitinfo$legend))
    splat(ledge)
  if(x$type != "rectangle")
    splat("Fraction of frame area:", signif(x$areafraction, 3))
  return(invisible(x))
}

as.data.frame.owin <- function(x, ..., drop=TRUE) {
  stopifnot(is.owin(x))
  switch(x$type,
         rectangle = { x <- as.polygonal(x) },
         polygonal = { },
         mask = {
           xy <- rasterxy.mask(x, drop=drop)
           if(!drop) xy <- append(xy, list(inside=as.vector(x$m)))
           return(as.data.frame(xy, ...))
         })
  b <- x$bdry
  ishole <- sapply(b, is.hole.xypolygon)
  sign <- (-1)^ishole
  b <- lapply(b, as.data.frame, ...)
  nb <- length(b)
  if(nb == 1)
    return(b[[1]])
  dfs <- mapply(cbind, b, id=as.list(seq_len(nb)), sign=as.list(sign),
                SIMPLIFY=FALSE)
  df <- do.call(rbind, dfs)
  return(df)
}

discretise <- function(X,eps=NULL,dimyx=NULL,xy=NULL) {
  verifyclass(X,"ppp")
  W <- X$window
  ok <- inside.owin(X$x,X$y,W)
  if(!all(ok))
    stop("There are points of X outside the window of X")
  all.null <- is.null(eps) & is.null(dimyx) & is.null(xy)
  if(W$type=="mask" & all.null) return(X)
  WM  <- as.mask(W,eps=eps,dimyx=dimyx,xy=xy)
  nok <- !inside.owin(X$x,X$y,WM)
  if(any(nok)) {
    ifix <- nearest.raster.point(X$x[nok],X$y[nok], WM)
    ifix <- cbind(ifix$row,ifix$col)
    WM$m[ifix] <- TRUE
  }
  X$window <- WM
  X
}

pixelcentres <- function (X, W=NULL,...) {
  X <- as.mask(as.owin(X), ...)
  if(is.null(W)) W <- as.rectangle(X)
  Y <- as.ppp(raster.xy(X,drop=TRUE),W=W)
  return(Y)
}

owin2polypath <- function(w) {
  w <- as.polygonal(w)
  b <- w$bdry
  xvectors <- lapply(b, getElement, name="x")
  yvectors <- lapply(b, getElement, name="y")
  xx <- unlist(lapply(xvectors, append, values=NA, after=FALSE))[-1]
  yy <- unlist(lapply(yvectors, append, values=NA, after=FALSE))[-1]
  return(list(x=xx, y=yy))
}

## generics which extract and assign the window of some object

Window <- function(X, ...) { UseMethod("Window") }

"Window<-" <- function(X, ..., value) { UseMethod("Window<-") }

