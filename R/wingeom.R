#
#	wingeom.S	Various geometrical computations in windows
#
#
#	$Revision: 4.82 $	$Date: 2013/11/01 06:49:45 $
#
#
#
#
#-------------------------------------

volume.owin <- function(x) { area.owin(x) }

area.owin <- function(w) {
	w <- as.owin(w)
        switch(w$type,
               rectangle = {
		width <- abs(diff(w$xrange))
		height <- abs(diff(w$yrange))
		area <- width * height
               },
               polygonal = {
                 area <- sum(unlist(lapply(w$bdry, area.xypolygon)))
               },
               mask = {
                 pixelarea <- abs(w$xstep * w$ystep)
                 npixels <- sum(w$m)
                 area <- pixelarea * npixels
               },
               stop("Unrecognised window type")
        )
        return(area)
}

perimeter <- function(w) {
  w <- as.owin(w)
  switch(w$type,
         rectangle = {
           return(2*(diff(w$xrange)+diff(w$yrange)))
         },
         polygonal={
           return(sum(lengths.psp(as.psp(w))))
         },
         mask={
           p <- as.polygonal(w)
           if(is.null(p)) return(NA)
           delta <- sqrt(w$xstep^2 + w$ystep^2)
           p <- simplify.owin(p, delta * 1.15)
           return(sum(lengths.psp(as.psp(p))))
         })
  return(NA)
}

sidelengths.owin <- function(x) {
  if(x$type != "rectangle")
    warning("Computing the side lengths of a non-rectangular window")
  with(x, c(diff(xrange), diff(yrange)))
}

shortside.owin <- function(x) { min(sidelengths(x)) }

eroded.areas <- function(w, r) {
  w <- as.owin(w)
  switch(w$type,
         rectangle = {
           width <- abs(diff(w$xrange))
           height <- abs(diff(w$yrange))
           areas <- pmax(width - 2 * r, 0) * pmax(height - 2 * r, 0)
         },
         polygonal = {
           # warning("Approximating polygonal window by digital image")
           w <- as.mask(w)
           areas <- eroded.areas(w, r)
         },
         mask = {
           # distances from each pixel to window boundary
           b <- bdist.pixels(w, style="matrix")
           # histogram breaks to satisfy hist()
           Bmax <- max(b, r)
           breaks <- c(-1,r,Bmax+1)
           # histogram of boundary distances
           h <- hist(b, breaks=breaks, plot=FALSE)$counts
           # reverse cumulative histogram
           H <- revcumsum(h)
           # drop first entry corresponding to r=-1
           H <- H[-1]
           # convert count to area
           pixarea <- w$xstep * w$ystep
           areas <- pixarea * H
         },
         stop("unrecognised window type")
         )
  areas
}	

even.breaks.owin <- function(w) {
	verifyclass(w, "owin")
        Rmax <- diameter(w)
        make.even.breaks(Rmax, Rmax/(100 * sqrt(2)))
}

unit.square <- function() { owin(c(0,1),c(0,1)) }

square <- function(r=1) {
  stopifnot(is.numeric(r))
  if(any(is.na(r) | !is.finite(r)))
    stop("argument r is NA or infinite")
  if(length(r) == 1) {
    stopifnot(r > 0)
    r <- c(0,r)
  } else if(length(r) == 2) {
    stopifnot(r[1] < r[2])
  } else stop("argument r must be a single number, or a vector of length 2")
  owin(r,r)
}

overlap.owin <- function(A, B) {
  # compute the area of overlap between two windows
  
  # check units
  if(!compatible.units(unitname(A), unitname(B)))
    warning("The two windows have incompatible units of length")
  
  At <- A$type
  Bt <- B$type
  if(At=="rectangle" && Bt=="rectangle") {
    xmin <- max(A$xrange[1],B$xrange[1])
    xmax <- min(A$xrange[2],B$xrange[2])
    if(xmax <= xmin) return(0)
    ymin <- max(A$yrange[1],B$yrange[1])
    ymax <- min(A$yrange[2],B$yrange[2])
    if(ymax <= ymin) return(0)
    return((xmax-xmin) * (ymax-ymin))
  }
  if((At=="rectangle" && Bt=="polygonal")
     || (At=="polygonal" && Bt=="rectangle")
     || (At=="polygonal" && Bt=="polygonal"))
  {
    AA <- as.polygonal(A)$bdry
    BB <- as.polygonal(B)$bdry
    area <- 0
    for(i in seq_along(AA))
      for(j in seq_along(BB))
        area <- area + overlap.xypolygon(AA[[i]], BB[[j]])
    return(area)
  }
  if(At=="mask") {
    # count pixels in A that belong to B
    pixelarea <- abs(A$xstep * A$ystep)
    x <- as.vector(raster.x(A)[A$m])
    y <- as.vector(raster.y(A)[A$m])
    ok <- inside.owin(x, y, B) 
    return(pixelarea * sum(ok))
  }
  if(Bt== "mask") {
    # count pixels in B that belong to A
    pixelarea <- abs(B$xstep * B$ystep)
    x <- as.vector(raster.x(B)[B$m])
    y <- as.vector(raster.y(B)[B$m])
    ok <- inside.owin(x, y, A)
    return(pixelarea * sum(ok))
  }
  stop("Internal error")
}

#
#  subset operator for window
#

"[.owin" <- 
function(x, i, ...) {
  if(!missing(i) && !is.null(i)) {
    if(is.im(i) && i$type == "logical") {
      # convert to window
      i <- as.owin(eval.im(ifelse1NA(i)))
    } else stopifnot(is.owin(i))
    x <- intersect.owin(x, i, fatal=FALSE)
  }
  return(x)
}
    
#
#
#  Intersection and union of windows
#
#

intersect.owin <- function(A, B, ..., fatal=TRUE) {
  liszt <- list(...)
  rasterinfo <- list()
  if(length(liszt) > 0) {
    # explicit arguments controlling raster info
    israster <- names(liszt) %in% names(formals(as.mask))
    rasterinfo <- liszt[israster]
    # handle intersection of more than two windows
    isowin <- unlist(lapply(liszt, is.owin))
    nextra <- sum(isowin)
    if(nextra > 0) {
      windows <- liszt[isowin]
      for(i in 1:nextra)
        B <- do.call("intersect.owin",
                     append(list(B, windows[[i]]), rasterinfo))
    }
  }
  if(missing(A) || is.null(A)) return(B)
  if(missing(B) || is.null(B)) return(A)
  verifyclass(A, "owin")
  verifyclass(B, "owin")
  #
  if(identical(A, B))
    return(A)

  # check units
  if(!compatible.units(unitname(A), unitname(B)))
    warning("The two windows have incompatible units of length")

  # determine intersection of x and y ranges
  xr <- intersect.ranges(A$xrange, B$xrange, fatal=fatal)
  yr <- intersect.ranges(A$yrange, B$yrange, fatal=fatal)
  if(!fatal && (is.null(xr) || is.null(yr)))
    return(NULL)
  C <- owin(xr, yr, unitname=unitname(A))

  if(is.empty(A) || is.empty(B))
    return(emptywindow(C))
           
  # Determine type of intersection
  
  Arect <- is.rectangle(A)
  Brect <- is.rectangle(B)
  Apoly <- is.polygonal(A)
  Bpoly <- is.polygonal(B)
  Amask <- is.mask(A)
  Bmask <- is.mask(B)
  
  # Rectangular case
  if(Arect && Brect)
    return(C)

  if(!Amask && !Bmask) {
    ####### Result is polygonal ############
    a <- lapply(as.polygonal(A)$bdry, reverse.xypolygon)
    b <- lapply(as.polygonal(B)$bdry, reverse.xypolygon)
    ab <- polyclip(a, b, "intersection", fillA="nonzero", fillB="nonzero")
    if(length(ab)==0)
      return(emptywindow(C))
    # ensure correct polarity
    totarea <- sum(unlist(lapply(ab, area.xypolygon)))
    if(totarea < 0)
      ab <- lapply(ab, reverse.xypolygon)
    AB <- owin(poly=ab, check=FALSE)
    AB <- rescue.rectangle(AB)
    return(AB)
  }

  ######### Result is a mask ##############
  
  # Restrict domain where possible 
  if(Arect)
    A <- C
  if(Brect)
    B <- C
  if(Amask)
    A <- trim.mask(A, C)
  if(Bmask)
    B <- trim.mask(B, C)

  # Did the user specify the pixel raster?
  if(length(rasterinfo) > 0) {
    # convert to masks with specified parameters, and intersect
    if(Amask) {
      A <- do.call("as.mask", append(list(A), rasterinfo))
      return(restrict.mask(A, B))
    } else {
      B <- do.call("as.mask", append(list(B), rasterinfo))
      return(restrict.mask(B, A))
    }
  } 
  
  # One mask and one rectangle?
  if(Arect && Bmask)
      return(B)
  if(Amask && Brect)
      return(A)

  # One mask and one polygon?
  if(Amask && !Bmask) 
    return(restrict.mask(A, B))
  if(!Amask && Bmask)
    return(restrict.mask(B, A))

  # Two existing masks?
  if(Amask && Bmask) {
    # choose the finer one
    if(A$xstep <= B$xstep)
      return(restrict.mask(A, B))
    else
      return(restrict.mask(B, A))
  }

  # No existing masks. No clipping applied so far.
  # Convert one window to a mask with default pixel raster, and intersect.
  if(Arect) {
    A <- as.mask(A)
    return(restrict.mask(A, B))
  } else {
    B <- as.mask(B)
    return(restrict.mask(B, A))
  }
}


union.owin <- function(A, B, ...) {
  liszt <- list(...)
  rasterinfo <- list()
  if(length(liszt) > 0) {
    # explicit arguments controlling raster info
    israster <- names(liszt) %in% names(formals(as.mask))
    rasterinfo <- liszt[israster]
    # handle intersection of more than two windows
    isowin <- unlist(lapply(liszt, is.owin))
    nextra <- sum(isowin)
    if(nextra > 0) {
      windows <- liszt[isowin]
      for(i in 1:nextra)
        B <- do.call("union.owin",
                     append(list(B, windows[[i]]), rasterinfo))
    }
  }
  #
  if(missing(A) || is.null(A) || is.empty(A)) return(B)
  if(missing(B) || is.null(B) || is.empty(B)) return(A)
  verifyclass(A, "owin")
  verifyclass(B, "owin")

  if(identical(A, B))
    return(A)

  # check units
  if(!compatible.units(unitname(A), unitname(B)))
    warning("The two windows have incompatible units of length")
  
  # Determine type of intersection
  
  Arect <- is.rectangle(A)
  Brect <- is.rectangle(B)
  Apoly <- is.polygonal(A)
  Bpoly <- is.polygonal(B)
  Amask <- is.mask(A)
  Bmask <- is.mask(B)

  # Result is not rectangular.
  # Create a rectangle to contain it.
  
  C <- owin(range(A$xrange, B$xrange),
            range(A$yrange, B$yrange),
            unitname=unitname(A))

  if(!Amask && !Bmask) {
    ####### Result is polygonal ############
    a <- lapply(as.polygonal(A)$bdry, reverse.xypolygon)
    b <- lapply(as.polygonal(B)$bdry, reverse.xypolygon)
    ab <- polyclip(a, b, "union", fillA="nonzero", fillB="nonzero")
    if(length(ab) == 0)
      return(emptywindow(C))
    # ensure correct polarity
    totarea <- sum(unlist(lapply(ab, area.xypolygon)))
    if(totarea < 0)
      ab <- lapply(ab, reverse.xypolygon)
    AB <- owin(poly=ab, check=FALSE)
    AB <- rescue.rectangle(AB)
    return(AB)
  }

  ####### Result is a mask ############

  # Determine pixel raster parameters
  if(length(rasterinfo) == 0) {
    rasterinfo <-
      if(Amask)
        list(xy=list(x=prolongseq(A$xcol, C$xrange),
               y=prolongseq(A$yrow, C$yrange)))
      else if(Bmask)
        list(xy=list(x=prolongseq(B$xcol, C$xrange),
               y=prolongseq(B$yrow, C$yrange)))
      else
        list()
  }

  # Convert C to mask
  C <- do.call("as.mask", append(list(w=C), rasterinfo))
      
  x <- as.vector(raster.x(C))
  y <- as.vector(raster.y(C))
  ok <- inside.owin(x, y, A) | inside.owin(x, y, B)

  if(all(ok)) {
    # result is a rectangle
    C <- as.rectangle(C)
  } else {
    # result is a mask
    C$m[] <- ok
  }
  return(C)
}

setminus.owin <- function(A, B, ...) {
  if(is.null(A) || is.empty(A)) return(B)
  if(is.null(B) || is.empty(B)) return(A)
  verifyclass(A, "owin")
  verifyclass(B, "owin")

  if(identical(A, B))
    return(emptywindow(as.rectangle(A)))

  # check units
  if(!compatible.units(unitname(A), unitname(B)))
    warning("The two windows have incompatible units of length")
  
  # Determine type of arguments
  
  Arect <- is.rectangle(A)
  Brect <- is.rectangle(B)
  Apoly <- is.polygonal(A)
  Bpoly <- is.polygonal(B)
  Amask <- is.mask(A)
  Bmask <- is.mask(B)

  # Case where A and B are both rectangular
  if(Arect && Brect) {
    C <- intersect.owin(A, B, fatal=FALSE)
    if(is.null(C)) return(A)
    return(complement.owin(C, A))
  }
    
  # Polygonal case

  if(!Amask && !Bmask) {
    ####### Result is polygonal ############
    a <- lapply(as.polygonal(A)$bdry, reverse.xypolygon)
    b <- lapply(as.polygonal(B)$bdry, reverse.xypolygon)
    ab <- polyclip(a, b, "minus", fillA="nonzero", fillB="nonzero")
    if(length(ab) == 0)
      return(emptywindow(C))
    # ensure correct polarity
    totarea <- sum(unlist(lapply(ab, area.xypolygon)))
    if(totarea < 0)
      ab <- lapply(ab, reverse.xypolygon)
    AB <- owin(poly=ab, check=FALSE)
    AB <- rescue.rectangle(AB)
    return(AB)
  }

  ####### Result is a mask ############

  # Determine pixel raster parameters
  rasterinfo <- 
    if((length(list(...)) > 0))
      list(...)
    else if(Amask)
      list(xy=list(x=A$xcol,
                   y=A$yrow))
    else if(Bmask)
      list(xy=list(x=B$xcol,
                   y=B$yrow))
    else
      list()

  # Convert A to mask
  AB <- do.call("as.mask", append(list(w=A), rasterinfo))
      
  x <- as.vector(raster.x(AB))
  y <- as.vector(raster.y(AB))
  ok <- inside.owin(x, y, A) & !inside.owin(x, y, B)

  if(!all(ok))
    AB$m[] <- ok
  else
    AB <- rescue.rectangle(AB)

  return(AB)
}

# auxiliary functions
  
trim.mask <- function(M, R, tolerant=TRUE) {
    # M is a mask,
    # R is a rectangle

    # Ensure R is a subset of bounding rectangle of M
    R <- owin(intersect.ranges(M$xrange, R$xrange),
              intersect.ranges(M$yrange, R$yrange))
    
    # Deal with very thin rectangles
    if(tolerant) {
      R$xrange <- adjustthinrange(R$xrange, M$xstep, M$xrange)
      R$yrange <- adjustthinrange(R$yrange, M$ystep, M$yrange)
    }

    # Extract subset of image grid
    within.range <- function(u, v) { (u >= v[1]) & (u <= v[2]) }
    yrowok <- within.range(M$yrow, R$yrange)
    xcolok <- within.range(M$xcol, R$xrange)
    if((ny <- sum(yrowok)) == 0 || (nx <- sum(xcolok)) == 0) 
      return(emptywindow(R))
    Z <- M
    Z$xrange <- R$xrange
    Z$yrange <- R$yrange
    Z$yrow <- M$yrow[yrowok]
    Z$xcol <- M$xcol[xcolok]
    Z$m <- M$m[yrowok, xcolok]
    if(ny < 2 || nx < 2)
      Z$m <- matrix(Z$m, nrow=ny, ncol=nx)
    Z$dim <- dim(Z$m)
    return(Z)
}

restrict.mask <- function(M, W) {
  # M is a mask, W is any window
  stopifnot(is.mask(M))
  stopifnot(inherits(W, "owin"))
  if(is.rectangle(W) == "rectangle")
    return(trim.mask(M, W))
  M <- trim.mask(M, as.rectangle(W))
  # Determine which pixels of M are inside W
  Mm <- M$m
  x <- as.vector(raster.x(M)[Mm])
  y <- as.vector(raster.y(M)[Mm])
  ok <- inside.owin(x, y, W)
  Mm[Mm] <- ok
  M$m <- Mm
  return(M)
}

# SUBSUMED IN rmhexpand.R
# expand.owin <- function(W, f=1) {
#
#  # expand bounding box of 'win'
#  # by factor 'f' in **area**
#  if(f <= 0)
#    stop("f must be > 0")
#  if(f == 1)
#    return(W)
#  bb <- bounding.box(W)
#  xr <- bb$xrange
#  yr <- bb$yrange
#  fff <- (sqrt(f) - 1)/2
#  Wexp <- owin(xr + fff * c(-1,1) * diff(xr),
#               yr + fff * c(-1,1) * diff(yr),
#               unitname=unitname(W))
#  return(Wexp)
#}

trim.rectangle <- function(W, xmargin=0, ymargin=xmargin) {
  if(!is.rectangle(W))
    stop("Internal error: tried to trim margin off non-rectangular window")
  xmargin <- ensure2vector(xmargin)
  ymargin <- ensure2vector(ymargin)
  if(any(xmargin < 0) || any(ymargin < 0))
    stop("values of xmargin, ymargin must be nonnegative")
  if(sum(xmargin) > diff(W$xrange))
    stop("window is too small to cut off margins of the width specified")
  if(sum(ymargin) > diff(W$yrange))
    stop("window is too small to cut off margins of the height specified")
  owin(W$xrange + c(1,-1) * xmargin,
       W$yrange + c(1,-1) * ymargin,
       unitname=unitname(W))
}

grow.rectangle <- function(W, xmargin=0, ymargin=xmargin) {
  xmargin <- ensure2vector(xmargin)
  ymargin <- ensure2vector(ymargin)
  if(any(xmargin < 0) || any(ymargin < 0))
    stop("values of xmargin, ymargin must be nonnegative")
  owin(W$xrange + c(-1,1) * xmargin,
       W$yrange + c(-1,1) * ymargin,
       unitname=unitname(W))
}


bdry.mask <- function(W) {
  verifyclass(W, "owin")
  W <- as.mask(W)
  m <- W$m
  nr <- nrow(m)
  nc <- ncol(m)
  b <-     (m != rbind(FALSE,       m[-nr, ]))
  b <- b | (m != rbind(m[-1, ], FALSE))
  b <- b | (m != cbind(FALSE,       m[, -nc]))
  b <- b | (m != cbind(m[, -1], FALSE))
  W$m <- b
  return(W)
}

vertices <- function(w) {
  verifyclass(w, "owin")
  if(is.empty(w))
    return(NULL)
  switch(w$type,
         rectangle={
           xr <- w$xrange
           yr <- w$yrange
           vert <- list(x=xr[c(1,2,2,1)], y=yr[c(1,1,2,2)])
         },
         polygonal={
           vert <- do.call("concatxy",w$bdry)
         },
         mask={
           b <- bdry.mask(w)
           xx <- raster.x(w)
           yy <- raster.y(w)
           vert <- list(x=as.vector(xx[b$m]),
                        y=as.vector(yy[b$m]))
         })
  return(vert)
}

diameter <- function(x) { UseMethod("diameter") }

diameter.owin <- function(x) {
  w <- as.owin(x)
  if(is.empty(w))
    return(NULL)
  vert <- vertices(w)
  if(length(vert$x) > 3) {
    # extract convex hull
    h <- with(vert, chull(x, y))
    vert <- with(vert, list(x=x[h], y=y[h]))
  }
  d <- pairdist(vert, squared=TRUE)
  return(sqrt(max(d)))
}

incircle <- function(W) {
  # computes the largest circle contained in W
  verifyclass(W, "owin")
  if(is.empty(W))
    return(NULL)
  if(is.rectangle(W)) {
    xr <- W$xrange
    yr <- W$yrange
    x0 <- mean(xr)
    y0 <- mean(yr)
    radius <- min(diff(xr), diff(yr))/2
    return(list(x=x0, y=y0, r=radius))
  }
  # compute distance to boundary
  D <- distmap(W, invert=TRUE)
  D <- D[W, drop=FALSE]
  # find maximum distance
  v <- D$v
  ok <- !is.na(v)
  Dvalues <- as.vector(v[ok])
  Dmax <- max(Dvalues)
  # find location of maximum
  locn <- which.max(Dvalues)
  locrow <- as.vector(row(v)[ok])[locn]
  loccol <- as.vector(col(v)[ok])[locn]
  x0 <- D$xcol[loccol]
  y0 <- D$yrow[locrow]
  if(is.mask(W)) {
    # radius could be one pixel diameter shorter than Dmax
    Dpixel <- sqrt(D$xstep^2 + D$ystep^2)
    radius <- max(0, Dmax - Dpixel)
  } else radius <- Dmax
  return(list(x=x0, y=y0, r=radius))
}

inpoint <- function(W) {
  # selects a point that is always inside the window.
  verifyclass(W, "owin")
  if(is.empty(W))
    return(NULL)
  if(is.rectangle(W))
    return(c(mean(W$xrange), mean(W$yrange)))
  if(is.polygonal(W)) {
    xy <- centroid.owin(W)
    if(inside.owin(xy$x, xy$y, W))
      return(xy)
  }
  W <- as.mask(W)
  Mm <- W$m
  Mrow <- as.vector(row(Mm)[Mm])
  Mcol <- as.vector(col(Mm)[Mm])
  selectmiddle <- function(x) { x[ceiling(length(x)/2)] }
  midcol <- selectmiddle(Mcol)
  midrow <- selectmiddle(Mrow[Mcol==midcol])
  x <- W$xcol[midcol]
  y <- W$yrow[midrow]
  return(c(x,y))
}

simplify.owin <- function(W, dmin) {
  verifyclass(W, "owin")
  if(is.empty(W))
    return(W)
  W <- as.polygonal(W)
  W$bdry <- lapply(W$bdry, simplify.xypolygon, dmin=dmin)
  return(W)
}

  
is.convex <- function(x) {
  verifyclass(x, "owin")
  if(is.empty(x))
    return(TRUE)
  switch(x$type,
         rectangle={return(TRUE)},
         polygonal={
           b <- x$bdry
           if(length(b) > 1)
             return(FALSE)
           b <- b[[1]]
           xx <- b$x
           yy <- b$y
           ch <- chull(xx,yy)
           return(length(ch) == length(xx))
         },
         mask={
           v <- vertices(x)
           v <- as.ppp(v, W=as.rectangle(x))
           ch <- convexhull.xy(v)
           edg <- as.psp(ch)
           edgedist <- nncross(v, edg, what="dist")
           pixdiam <- sqrt(x$xstep^2 + x$ystep^2)
           return(all(edgedist <= pixdiam))
         })
  return(as.logical(NA))
}

convexhull <- function(x) {
  if(inherits(x, "owin")) 
    v <- vertices(x)
  else if(inherits(x, "psp"))
    v <- endpoints.psp
  else if(inherits(x, "ppp"))
    v <- x
  else {
    x <- as.owin(x)
    v <- vertices(x)
  }
  b <- as.rectangle(x)
  if(is.empty(x))
    return(emptywindow(b))
  ch <- convexhull.xy(v)
  out <- rebound.owin(ch, b)
  return(out)
}

  
is.empty <- function(x) { UseMethod("is.empty") }

is.empty.default <- function(x) { length(x) == 0 }

is.empty.owin <- function(x) {
  switch(x$type,
         rectangle=return(FALSE),
         polygonal=return(length(x$bdry) == 0),
         mask=return(!any(x$m)))
  return(NA)
}

emptywindow <- function(w) {
  w <- as.owin(w)
  out <- owin(w$xrange, w$yrange, poly=list(), unitname=unitname(w))
  return(out)
}

