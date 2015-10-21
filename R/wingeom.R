#
#	wingeom.S	Various geometrical computations in windows
#
#
#	$Revision: 4.104 $	$Date: 2015/10/21 09:06:57 $
#
#
#
#
#-------------------------------------

volume.owin <- function(x) { area.owin(x) }

area <- function(w) UseMethod("area")

area.default <- function(w) area.owin(as.owin(w))

area.owin <- function(w) {
    stopifnot(is.owin(w))
        switch(w$type,
               rectangle = {
		width <- abs(diff(w$xrange))
		height <- abs(diff(w$yrange))
		area <- width * height
               },
               polygonal = {
                 area <- sum(unlist(lapply(w$bdry, Area.xypolygon)))
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
           return(sum(lengths.psp(edges(w))))
         },
         mask={
           p <- as.polygonal(w)
           if(is.null(p)) return(NA)
           delta <- sqrt(w$xstep^2 + w$ystep^2)
           p <- simplify.owin(p, delta * 1.15)
           return(sum(lengths.psp(edges(p))))
         })
  return(NA)
}

framebottomleft <- function(w) {
  f <- Frame(w)
  c(f$xrange[1], f$yrange[1])
}

sidelengths.owin <- function(x) {
  if(x$type != "rectangle")
    warning("Computing the side lengths of a non-rectangular window")
  with(x, c(diff(xrange), diff(yrange)))
}

shortside.owin <- function(x) { min(sidelengths(x)) }

eroded.areas <- function(w, r, subset=NULL) {
  w <- as.owin(w)
  if(!is.null(subset) && !is.mask(w))
    w <- as.mask(w)
  switch(w$type,
         rectangle = {
           width <- abs(diff(w$xrange))
           height <- abs(diff(w$yrange))
           areas <- pmax(width - 2 * r, 0) * pmax(height - 2 * r, 0)
         },
         polygonal = {
           ## warning("Approximating polygonal window by digital image")
           w <- as.mask(w)
           areas <- eroded.areas(w, r)
         },
         mask = {
           ## distances from each pixel to window boundary
           b <- if(is.null(subset)) bdist.pixels(w, style="matrix") else 
                bdist.pixels(w)[subset, drop=TRUE, rescue=FALSE]
           ## histogram breaks to satisfy hist()
           Bmax <- max(b, r)
           breaks <- c(-1,r,Bmax+1)
           ## histogram of boundary distances
           h <- hist(b, breaks=breaks, plot=FALSE)$counts
           ## reverse cumulative histogram
           H <- revcumsum(h)
           ## drop first entry corresponding to r=-1
           H <- H[-1]
           ## convert count to area
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

square <- function(r=1, unitname=NULL) {
  stopifnot(is.numeric(r))
  if(any(is.na(r) | !is.finite(r)))
    stop("argument r is NA or infinite")
  if(length(r) == 1) {
    stopifnot(r > 0)
    r <- c(0,r)
  } else if(length(r) == 2) {
    stopifnot(r[1] < r[2])
  } else stop("argument r must be a single number, or a vector of length 2")
  owin(r,r, unitname=unitname)
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
    rxy <- rasterxy.mask(A, drop=TRUE)
    x <- rxy$x
    y <- rxy$y
    ok <- inside.owin(x, y, B) 
    return(pixelarea * sum(ok))
  }
  if(Bt== "mask") {
    # count pixels in B that belong to A
    pixelarea <- abs(B$xstep * B$ystep)
    rxy <- rasterxy.mask(B, drop=TRUE)
    x <- rxy$x
    y <- rxy$y
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

intersect.owin <- function(..., fatal=TRUE, p) {
  argh <- list(...)
  ## p is a list of arguments to polyclip::polyclip
  if(missing(p) || is.null(p)) p <- list()
  ## handle 'solist' objects
  argh <- expandSpecialLists(argh, "solist")
  rasterinfo <- list()
  if(length(argh) > 0) {
    # explicit arguments controlling raster info
    israster <- names(argh) %in% names(formals(as.mask))
    if(any(israster)) {
      rasterinfo <- argh[israster]
      ## remaining arguments
      argh <- argh[!israster]
    }
  }
  ## look for window arguments
  isowin <- sapply(argh, is.owin)
  if(any(!isowin))
    warning("Some arguments were not windows")
  argh <- argh[isowin]
  nwin <- length(argh)

  if(nwin == 0) {
    warning("No windows were given")
    return(NULL)
  }

  ## at least one window
  A <- argh[[1]]
  if(nwin == 1) return(A)
  ## at least two windows
  B <- argh[[2]]
  
  if(nwin > 2) {
    ## handle union of more than two windows
    windows <- argh[-c(1,2)]
    ## determine a common set of parameters for polyclip
    p <- commonPolyclipArgs(A, B, do.call(boundingbox, windows), p=p)
    ## absorb all windows into B
    for(i in seq_along(windows))
      B <- do.call(intersect.owin,
                   append(list(B, windows[[i]], p=p),
                            rasterinfo))
  }

  ## There are now only two windows
  if(is.empty(A)) return(A)
  if(is.empty(B)) return(B)

  if(identical(A, B))
    return(A)

  # check units
  if(!compatible(unitname(A), unitname(B)))
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
#  Apoly <- is.polygonal(A)
#  Bpoly <- is.polygonal(B)
  Amask <- is.mask(A)
  Bmask <- is.mask(B)
  
  # Rectangular case
  if(Arect && Brect)
    return(C)

  if(!Amask && !Bmask) {
    ####### Result is polygonal ############
    a <- lapply(as.polygonal(A)$bdry, reverse.xypolygon)
    b <- lapply(as.polygonal(B)$bdry, reverse.xypolygon)
    ab <- do.call(polyclip::polyclip,
                  append(list(a, b, "intersection",
                              fillA="nonzero", fillB="nonzero"),
                         p))
    if(length(ab)==0)
      return(emptywindow(C))
    # ensure correct polarity
    totarea <- sum(unlist(lapply(ab, Area.xypolygon)))
    if(totarea < 0)
      ab <- lapply(ab, reverse.xypolygon)
    AB <- owin(poly=ab, check=FALSE, unitname=unitname(A))
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


union.owin <- function(..., p) {
  argh <- list(...)
  ## weed out NULL arguments
  argh <- argh[!sapply(argh, is.null)]
  ## p is a list of arguments to polyclip::polyclip
  if(missing(p) || is.null(p)) p <- list()
  ## handle 'solist' objects
  argh <- expandSpecialLists(argh, "solist")
  rasterinfo <- list()
  if(length(argh) > 0) {
    ## arguments controlling raster info
    israster <- names(argh) %in% names(formals(as.mask))
    if(any(israster)) {
      rasterinfo <- argh[israster]
      ## remaining arguments
      argh <- argh[!israster]
    }
  }
  ## look for window arguments
  isowin <- sapply(argh, is.owin)
  if(any(!isowin))
    warning("Some arguments were not windows")
  argh <- argh[isowin]
  ## 
  nwin <- length(argh)
  if(nwin == 0) {
    warning("No windows were given")
    return(NULL)
  }
  ## find non-empty ones
  if(any(isemp <- sapply(argh, is.empty)))
    argh <- argh[!isemp]
  nwin <- length(argh)
  if(nwin == 0) {
    warning("All windows were empty")
    return(NULL)
  }
  ## at least one window
  A <- argh[[1]]
  if(nwin == 1) return(A)

  ## more than two windows
  if(nwin > 2) {
    ## check if we need polyclip
    somepoly <- !all(sapply(argh, is.mask))
    if(somepoly) {
      ## determine a common set of parameters for polyclip
      p <- commonPolyclipArgs(do.call(boundingbox, argh), p=p)
      ## apply these parameters now to avoid numerical errors
      argh <- applyPolyclipArgs(argh, p=p)
      A <- argh[[1]]
    } 
    ## absorb all windows into A without rescaling
    nullp <- list(eps=1, x0=0, y0=0)
    for(i in 2:nwin) 
      A <- do.call(union.owin,
                   append(list(A, argh[[i]], p=nullp),
                          rasterinfo))
    if(somepoly) {
      ## undo rescaling
      A <- reversePolyclipArgs(A, p=p)
    }
    return(A)
  }

  ## Exactly two windows
  B <- argh[[2]]
  if(identical(A, B))
    return(A)

  ## check units
  if(!compatible(unitname(A), unitname(B)))
    warning("The two windows have incompatible units of length")
  
  ## Determine type of intersection
  
##  Arect <- is.rectangle(A)
##  Brect <- is.rectangle(B)
##  Apoly <- is.polygonal(A)
##  Bpoly <- is.polygonal(B)
  Amask <- is.mask(A)
  Bmask <- is.mask(B)

  ## Create a rectangle to contain the result
  
  C <- owin(range(A$xrange, B$xrange),
            range(A$yrange, B$yrange),
            unitname=unitname(A))

  if(!Amask && !Bmask) {
    ####### Result is polygonal (or rectangular) ############
    a <- lapply(as.polygonal(A)$bdry, reverse.xypolygon)
    b <- lapply(as.polygonal(B)$bdry, reverse.xypolygon)
    ab <- do.call(polyclip::polyclip,
                  append(list(a, b, "union",
                              fillA="nonzero", fillB="nonzero"),
                         p))
    if(length(ab) == 0)
      return(emptywindow(C))
    ## ensure correct polarity
    totarea <- sum(unlist(lapply(ab, Area.xypolygon)))
    if(totarea < 0)
      ab <- lapply(ab, reverse.xypolygon)
    AB <- owin(poly=ab, check=FALSE, unitname=unitname(A))
    AB <- rescue.rectangle(AB)
    return(AB)
  }

  ####### Result is a mask ############

  ## Determine pixel raster parameters
  if(length(rasterinfo) == 0) {
    rasterinfo <-
      if(Amask)
        list(xy=list(x=as.numeric(prolongseq(A$xcol, C$xrange)),
               y=as.numeric(prolongseq(A$yrow, C$yrange))))
      else if(Bmask)
        list(xy=list(x=as.numeric(prolongseq(B$xcol, C$xrange)),
               y=as.numeric(prolongseq(B$yrow, C$yrange))))
      else
        list()
  }

  ## Convert C to mask
  C <- do.call("as.mask", append(list(w=C), rasterinfo))

  rxy <- rasterxy.mask(C)
  x <- rxy$x
  y <- rxy$y
  ok <- inside.owin(x, y, A) | inside.owin(x, y, B)

  if(all(ok)) {
    ## result is a rectangle
    C <- as.rectangle(C)
  } else {
    ## result is a mask
    C$m[] <- ok
  }
  return(C)
}

setminus.owin <- function(A, B, ..., p) {
  if(is.null(A) || is.empty(A)) return(B)
  if(is.null(B) || is.empty(B)) return(A)
  verifyclass(A, "owin")
  verifyclass(B, "owin")
  if(identical(A, B))
    return(emptywindow(as.rectangle(A)))

  ## p is a list of arguments to polyclip::polyclip
  if(missing(p) || is.null(p)) p <- list()

  ## check units
  if(!compatible(unitname(A), unitname(B)))
    warning("The two windows have incompatible units of length")
  
  ## Determine type of arguments
  
  Arect <- is.rectangle(A)
  Brect <- is.rectangle(B)
##  Apoly <- is.polygonal(A)
##  Bpoly <- is.polygonal(B)
  Amask <- is.mask(A)
  Bmask <- is.mask(B)

  ## Case where A and B are both rectangular
  if(Arect && Brect) {
    C <- intersect.owin(A, B, fatal=FALSE)
    if(is.null(C)) return(A)
    return(complement.owin(C, A))
  }
    
  ## Polygonal case

  if(!Amask && !Bmask) {
    ####### Result is polygonal ############
    a <- lapply(as.polygonal(A)$bdry, reverse.xypolygon)
    b <- lapply(as.polygonal(B)$bdry, reverse.xypolygon)
    ab <- do.call(polyclip::polyclip,
                  append(list(a, b, "minus",
                              fillA="nonzero", fillB="nonzero"),
                         p))
    if(length(ab) == 0)
      return(emptywindow(A))
    ## ensure correct polarity
    totarea <- sum(unlist(lapply(ab, Area.xypolygon)))
    if(totarea < 0)
      ab <- lapply(ab, reverse.xypolygon)
    AB <- owin(poly=ab, check=FALSE, unitname=unitname(A))
    AB <- rescue.rectangle(AB)
    return(AB)
  }

  ####### Result is a mask ############

  ## Determine pixel raster parameters
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

  ## Convert A to mask
  AB <- do.call("as.mask", append(list(w=A), rasterinfo))

  rxy <- rasterxy.mask(AB)
  x <- rxy$x
  y <- rxy$y
  ok <- inside.owin(x, y, A) & !inside.owin(x, y, B)

  if(!all(ok))
    AB$m[] <- ok
  else
    AB <- rescue.rectangle(AB)

  return(AB)
}

## auxiliary functions

commonPolyclipArgs <- function(..., p=NULL) {
  # compute a common resolution for polyclip operations
  # on several windows
  if(!is.null(p) && !is.null(p$eps) && !is.null(p$x0) && !is.null(p$y0))
    return(p)
  bb <- boundingbox(...)
  xr <- bb$xrange
  yr <- bb$yrange
  eps <- p$eps %orifnull% max(diff(xr), diff(yr))/(2^31)
  x0  <- p$x0  %orifnull% mean(xr)
  y0  <- p$y0  %orifnull% mean(yr)
  return(list(eps=eps, x0=x0, y0=y0))
}

applyPolyclipArgs <- function(x, p=NULL) {
  if(is.null(p)) return(x)
  y <- lapply(x, shift, vec=-c(p$x0, p$y0))
  z <- lapply(y, scalardilate, f=1/p$eps)
  return(z)
}

reversePolyclipArgs <- function(x, p=NULL) {
  if(is.null(p)) return(x)
  y <- scalardilate(x, f=p$eps)
  z <- shift(y, vec=c(p$x0, p$y0))
  return(z)
}

trim.mask <- function(M, R, tolerant=TRUE) {
    ## M is a mask,
    ## R is a rectangle

    ## Ensure R is a subset of bounding rectangle of M
    R <- owin(intersect.ranges(M$xrange, R$xrange),
              intersect.ranges(M$yrange, R$yrange))
    
    ## Deal with very thin rectangles
    if(tolerant) {
      R$xrange <- adjustthinrange(R$xrange, M$xstep, M$xrange)
      R$yrange <- adjustthinrange(R$yrange, M$ystep, M$yrange)
    }

    ## Extract subset of image grid
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
  ## M is a mask, W is any window
  stopifnot(is.mask(M))
  stopifnot(inherits(W, "owin"))
  if(is.rectangle(W))
    return(trim.mask(M, W))
  M <- trim.mask(M, as.rectangle(W))
  ## Determine which pixels of M are inside W
  rxy <- rasterxy.mask(M, drop=TRUE)
  x <- rxy$x
  y <- rxy$y
  ok <- inside.owin(x, y, W)
  Mm <- M$m
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
#  bb <- boundingbox(W)
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

grow.mask <- function(M, xmargin=0, ymargin=xmargin) {
  stopifnot(is.mask(M))
  m <- as.matrix(M)
  Rplus <- grow.rectangle(as.rectangle(M), xmargin, ymargin)
  ## extend the raster
  xcolplus <- prolongseq(M$xcol, Rplus$xrange)
  yrowplus <- prolongseq(M$yrow, Rplus$yrange)
  mplus <- matrix(FALSE, length(yrowplus), length(xcolplus))
  ## pad out the mask entries
  nleft <- attr(xcolplus, "nleft")
  nright <- attr(xcolplus, "nright")
  nbot <- attr(yrowplus, "nleft")
  ntop <- attr(yrowplus, "nright")
  mplus[ (nleft+1):(length(yrowplus)-nright),
         (nbot+1):(length(xcolplus)-ntop) ] <- m
  ## pack up
  result <- owin(xrange=Rplus$xrange,
                 yrange=Rplus$yrange,
                 xcol=as.numeric(xcolplus),
                 yrow=as.numeric(yrowplus),
                 mask=mplus,
                 unitname=unitname(M))
  return(result)
}
  
bdry.mask <- function(W) {
  verifyclass(W, "owin")
  W <- as.mask(W)
  m <- W$m
  nr <- nrow(m)
  nc <- ncol(m)
  if(!spatstat.options('Cbdrymask')) {
    ## old interpreted code
    b <-     (m != rbind(FALSE,       m[-nr, ]))
    b <- b | (m != rbind(m[-1, ], FALSE))
    b <- b | (m != cbind(FALSE,       m[, -nc]))
    b <- b | (m != cbind(m[, -1], FALSE))
  } else {
    b <- integer(nr * nc)
    z <- .C("bdrymask",
            nx = as.integer(nc),
            ny = as.integer(nr),
            m = as.integer(m),
            b = as.integer(b))
    b <- matrix(as.logical(z$b), nr, nc)
  }
  W$m <- b
  return(W)
}

vertices <- function(w) {
  UseMethod("vertices")
}

vertices.owin <- function(w) {
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
           bm <- bdry.mask(w)$m
           rxy <- rasterxy.mask(w)
           xx <- rxy$x
           yy <- rxy$y
           vert <- list(x=as.vector(xx[bm]),
                        y=as.vector(yy[bm]))
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

##    radius of inscribed circle

inradius <- function(W) {
  stopifnot(is.owin(W))
  if(W$type == "rectangle") diameter(W)/2 else max(distmap(W, invert=TRUE))
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
           edg <- edges(ch)
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
    v <- endpoints.psp(x)
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

discs <- function(centres, radii=marks(centres)/2, ..., 
                  separate=FALSE, mask=FALSE, trim=TRUE, delta=NULL) {
  stopifnot(is.ppp(centres))
  n <- npoints(centres)
  if(n == 0) return(emptywindow(Frame(centres)))
  check.nvector(radii, npoints(centres))
  stopifnot(all(radii > 0))
  if(is.null(delta)) delta <- 2 * pi * min(radii)/16
  if(separate) {
    D <- list()
    W <- disc(centre=centres[1], radius=radii[1], delta=delta)
    D[[1]] <- W
    if(n == 1) return(D)
    for(i in 2:n) 
      D[[i]] <- disc(centre=centres[i], radius=radii[i], delta=delta)
    return(D)
  } else if(mask) {
    M <- as.mask(Window(centres), ...)
    z <- .C("discs2grid",
            nx    = as.integer(M$dim[2]),
            x0    = as.double(M$xcol[1]),
            xstep = as.double(M$xstep),  
            ny    = as.integer(M$dim[1]),
            y0    = as.double(M$yrow[1]),
            ystep = as.double(M$ystep), 
            nd    = as.integer(n),
            xd    = as.double(centres$x),
            yd    = as.double(centres$y),
            rd    = as.double(radii), 
            out   = as.integer(integer(prod(M$dim))))
    M$m[] <- as.logical(z$out)
    return(M)
  } else {
    W <- disc(centre=centres[1], radius=radii[1], delta=delta)
    if(n == 1) return(W)
    for(i in 2:n) {
      Di <- disc(centre=centres[i], radius=radii[i], delta=delta)
      W <- union.owin(W, Di)
    }
    if(trim) W <- intersect.owin(W, Window(centres))
    return(W)
  }
}

harmonise.owin <- harmonize.owin <- function(...) {
  argz <- list(...)
  wins <- solapply(argz, as.owin)
  if(length(wins) < 2L) return(wins)
  ismask <- sapply(wins, is.mask)
  if(!any(ismask)) return(wins)
  comgrid <- do.call(commonGrid, lapply(argz, as.owin))
  result <- solapply(argz, "[", i=comgrid, drop=FALSE)
  return(result)
}

