#
#	dummy.S
#
#	Utilities for generating patterns of dummy points
#
#       $Revision: 5.31 $     $Date: 2015/10/21 09:06:57 $
#
#	corners()	corners of window
#	gridcenters()	points of a rectangular grid
#	stratrand()	random points in each tile of a rectangular grid
#	spokes()	Rolf's 'spokes' arrangement
#	
#	concatxy()	concatenate any lists of x, y coordinates
#
#	default.dummy()	Default action to create a dummy pattern
#		
	
corners <- function(window) {
	window <- as.owin(window)
	x <- window$xrange[c(1,2,1,2)]
	y <- window$yrange[c(1,1,2,2)]
	return(list(x=x, y=y))
}

gridcenters <-	
gridcentres <- function(window, nx, ny) {
	window <- as.owin(window)
	xr <- window$xrange
	yr <- window$yrange
	x <- seq(from=xr[1], to=xr[2], length.out = 2 * nx + 1)[2 * (1:nx)]
	y <- seq(from=yr[1], to=yr[2], length.out = 2 * ny + 1)[2 * (1:ny)]
	x <- rep.int(x, ny)
	y <- rep.int(y, rep.int(nx, ny))
	return(list(x=x, y=y))
}

stratrand <- function(window,nx,ny, k=1) {
	
	# divide window into an nx * ny grid of tiles
	# and place k points at random in each tile
	
	window <- as.owin(window)

	wide  <- diff(window$xrange)/nx
	high  <- diff(window$yrange)/ny
        cent <- gridcentres(window, nx, ny)
	cx <- rep.int(cent$x, k)
	cy <- rep.int(cent$y, k)
	n <- nx * ny * k
	x <- cx + runif(n, min = -wide/2, max = wide/2)
	y <- cy + runif(n, min = -high/2, max = high/2)
	return(list(x=x,y=y))
}

tilecentroids <- function (W, nx, ny)
{
  W <- as.owin(W)
  if(W$type == "rectangle")
    return(gridcentres(W, nx, ny))
  else {
    # approximate
    W   <- as.mask(W)
    rxy <- rasterxy.mask(W, drop=TRUE)
    xx  <- rxy$x
    yy  <- rxy$y
    pid <- gridindex(xx,yy,W$xrange,W$yrange,nx,nx)$index
    x   <- tapply(xx,pid,mean)
    y   <- tapply(yy,pid,mean)
    return(list(x=x,y=y))
  }
}

cellmiddles <- local({
  # auxiliary 
  middle <- function(v) { n <- length(v);
                          mid <- ceiling(n/2);
                          v[mid]}

  dcut <- function(x, nx, xrange) {
    dx <- diff(xrange)/nx
    fx <- ((x - xrange[1])/dx) %% 1
    bx <- dx * pmin(fx, 1-fx)
    bx
  }
  
  # main
  cellmiddles <- function (W, nx, ny, npix=NULL, distances=FALSE) {
    if(W$type == "rectangle")
      return(gridcentres(W, nx, ny))

    # pixel approximation to window
    # This matches the pixel approximation used to compute tile areas
    # and ensures that dummy points are generated only inside those tiles
    # that have nonzero digital area
    M   <- as.mask(W, dimyx=rev(npix))
    xx <- as.vector(rasterx.mask(M, drop=TRUE))
    yy <- as.vector(rastery.mask(M, drop=TRUE))
    pid <- gridindex(xx,yy,W$xrange,W$yrange,nx,ny)$index

    # compute tile centroids
    xmid <- tapply(xx, pid, mean)
    ymid <- tapply(yy, pid, mean)
    # check whether they are inside window
    ok <- inside.owin(xmid, ymid, W)
    if(all(ok))
      return(list(x=xmid, y=ymid))

    # some problem tiles
    bad <- rep.int(TRUE, nx * ny)
    bad[as.integer(names(xmid))] <- !ok
    badpid <- bad[pid]
    if(!distances) {
       midpix <- tapply(seq_along(pid)[badpid], pid[badpid], middle)
    } else {
      # find 'middle' points using boundary distances
      Dlines <- im(outer(dcut(M$yrow,ny,M$yrange),
                         dcut(M$xcol,nx,M$xrange),
                         "pmin"),
                   M$xcol, M$yrow, M$xrange, M$yrange)
      Dbdry <- bdist.pixels(M)
      Dtile <- eval.im(pmin(Dlines, Dbdry))
      dtile <- as.vector(Dtile[M])
      df <- data.frame(dtile=dtile, id=seq_along(dtile))[badpid, , drop=FALSE]
      midpix <- by(df, pid[badpid], midpixid)
    }
    xmid[!ok] <- xx[midpix]
    ymid[!ok] <- yy[midpix]
    return(list(x=xmid,y=ymid))
  }

  midpixid <- function(z) { z$id[which.max(z$dtile)] }
  
  cellmiddles
})

spokes <- function(x, y, nrad = 3, nper = 3, fctr = 1.5, Mdefault=1) {
	#
	# Rolf Turner's "spokes" arrangement
	#
	# Places dummy points on radii of circles 
	# emanating from each data point x[i], y[i]
	#
	#       nrad:    number of radii from each data point
	#       nper:	 number of dummy points per radius
	#       fctr:	 length of largest radius = fctr * M
	#                where M is mean nearest neighbour distance in data
	#
        pat <- inherits(x,"ppp")
        if(pat) w <- x$w
        if(checkfields(x,c("x","y"))) {
          y <- x$y
          x <- x$x
        }
        M <- if(length(x) > 1) mean(nndist(x,y)) else Mdefault
	lrad  <- fctr * M / nper
	theta <- 2 * pi * (1:nrad)/nrad
	cs    <- cos(theta)
	sn    <- sin(theta)
	xt    <- lrad * as.vector((1:nper) %o% cs)
	yt    <- lrad * as.vector((1:nper) %o% sn)
	xd    <- as.vector(outer(x, xt, "+"))
	yd    <- as.vector(outer(y, yt, "+"))
	
        tmp <- list(x = xd, y = yd)
        if(pat) return(as.ppp(tmp,W=w)[w]) else return(tmp)
}
	
# concatenate any number of list(x,y) into a list(x,y)
		
concatxy <- function(...) {
	x <- unlist(lapply(list(...), getElement, name="x"))
	y <- unlist(lapply(list(...), getElement, name="y"))
	if(length(x) != length(y))
		stop("Internal error: lengths of x and y unequal")
	return(list(x=x,y=y))
}

#------------------------------------------------------------

default.dummy <- function(X, nd=NULL, random=FALSE, ntile=NULL, npix = NULL,
                          quasi=FALSE, ..., eps=NULL, verbose=FALSE) {
  # default action to create dummy points.
  # regular grid of nd[1] * nd[2] points
  # plus corner points of window frame,
  # all clipped to window.
  X <- as.ppp(X)
  win <- X$window
  #
  # default dimensions
  a <- default.n.tiling(X, nd=nd, ntile=ntile, npix=npix,
                        eps=eps, random=random, quasi=quasi, verbose=verbose)
  nd    <- a$nd
  ntile <- a$ntile
  npix  <- a$npix
  periodsample <- !quasi && !random &&
                  is.mask(win) &&
                  all(nd %% win$dim == 0)
  # make dummy points
  dummy <- if(quasi) rQuasi(prod(nd), as.rectangle(win)) else
           if(random) stratrand(win, nd[1], nd[2], 1) else 
           cellmiddles(win, nd[1], nd[2], npix)
  dummy <- as.ppp(dummy, win, check=FALSE)
  # restrict to window
  if(!is.rectangle(win) && !periodsample)
    dummy <- dummy[win]
  # corner points
  corn <- as.ppp(corners(win), win, check=FALSE)
  corn <- corn[win]
  dummy <- superimpose(dummy, corn, W=win, check=FALSE)
  if(dummy$n == 0)
    stop("None of the dummy points lies inside the window")
  # pass parameters for computing weights
  attr(dummy, "dummy.parameters") <-
    list(nd=nd, random=random, quasi=quasi, verbose=verbose)
  attr(dummy, "weight.parameters") <-
    append(list(...), list(ntile=ntile, verbose=verbose, npix=npix))
  return(dummy)
}


# Criteria:
#   for rectangular windows,
#       R1.  nd >= ntile
#   for non-rectangular windows,
#       R2. nd should be a multiple of ntile
#       R3. each dummy point is also a pixel of the npix grid
#       R4. npix should ideally be a multiple of nd, for speed
#       R5. npix should be large, for accuracy
#       R6. npix should not be too large, for speed
#       R7. if the window is a mask, npix should ideally be
#           a multiple of the mask array dimensions, for speed.
#

default.n.tiling <- local({
  # auxiliary
  ensure2print <- function(x, verbose=TRUE, blah="user specified") {
    xname <- short.deparse(substitute(x))
    x <- ensure2vector(x)
    if(verbose)
      cat(paste(blah, xname, "=", x[1], "*", x[2], "\n"))
    x
  }
  minmultiple <- function(n, lo, hi) {
    if(lo > hi) {
      temp <- hi
      hi <- lo
      lo <- temp
    }
    if(n > hi) return(hi)
    m <- n * (floor(lo/n):ceiling(hi/n))
    m <- m[m >= n & m >= lo & m <= hi]
    if(length(m) > 0) min(m) else hi
  }
    
  mindivisor <- function(N, lo, Nbig) {
    d <- divisors(N)
    ok <- (d >= lo)
    if(any(ok)) return(min(d[ok]))
    m <- floor(Nbig/N)
    d <- unlist(lapply(as.list(seq_len(m) * N), divisors))
    d <- sort(unique(d))
    ok <- (d >= lo)
    if(any(ok)) return(min(d[ok]))
    return(Nbig)
  }

  min2mul <- function(n, lo, hi) 
    c(minmultiple(n[1], lo[1], hi[1]),
      minmultiple(n[2], lo[2], hi[2]))

  min2div <- function(N, lo, Nbig) 
    c(mindivisor(N[1], lo[1], Nbig[1]),
      mindivisor(N[2], lo[2], Nbig[2]))

  maxdiv <- function(n, k=1) {
    if(length(n) > 1)
      return(c(maxdiv(n[1], k),
               maxdiv(n[2], k)))
    ## k-th largest divisor other than n
    d <- divisors(n)
    m <- length(d)
    ans <- if(m == 2) n else if(m < 2+k) d[2] else d[m-k]
    return(ans)
  }

  # main
  default.n.tiling <- function(X,
                               nd=NULL, ntile=NULL, npix=NULL,
                               eps=NULL,
                               random=FALSE, quasi=FALSE, verbose=TRUE) {
  # computes dimensions of rectangular grids of 
  #     - dummy points  (nd) (eps)
  #     - tiles for grid weights (ntile)
  #     - pixels for approximating area (npix)
  # for data pattern X.
  #
  verifyclass(X, "ppp")
  win <- X$window
  pixels <- (win$type != "rectangle")
  
  if(nd.given <- !is.null(nd)) 
    nd <- ensure2print(nd, verbose)
  if(ntile.given <- !is.null(ntile)) 
    ntile <- ensure2print(ntile, verbose)
  if(npix.given <- !is.null(npix)) 
    npix <- ensure2print(npix, verbose)

  if(pixels) 
    sonpixel <- rev(ensure2print(spatstat.options("npixel"), verbose, ""))

  ndummy.min <- ensure2print(spatstat.options("ndummy.min"), verbose, "")
  ndminX <- pmax(ndummy.min, 10 * ceiling(2 * sqrt(X$n)/10))
  ndminX <- ensure2vector(ndminX)

  if(!is.null(eps)) {
    eps <- ensure2print(eps, verbose)
    Xbox <- as.rectangle(as.owin(X))
    sides <- with(Xbox, c(diff(xrange), diff(yrange)))
    ndminX <- pmax(ndminX, ceiling(sides/eps))
  }

  # range of acceptable values for npix
  if(npix.given)
    Nmin <- Nmax <- npix
  else 
    switch(win$type,
           rectangle = {
             Nmin <- ensure2vector(X$n)
             Nmax <- Inf
           },
           polygonal = {
             Nmin <- sonpixel
             Nmax <- 4 * sonpixel
           },
           mask={
             nmask <- rev(win$dim)
             Nmin <- nmask
             Nmax <- pmax(2 * nmask, 4 * sonpixel)
           })

  # determine values of nd and ntile

  if(nd.given && !ntile.given) {
    # ntile must be a divisor of nd
    if(any(nd > Nmax))
      warning("number of dummy points nd exceeds maximum pixel dimensions")
    ntile <- min2div(nd, ndminX, nd)
  } else if(!nd.given && ntile.given) {
    # nd must be a multiple of ntile
    nd <- min2mul(ntile, ndminX, Nmin)
    if(any(nd >= Nmin))
      nd <- ntile
  } else if(!nd.given && !ntile.given) {
     if(!pixels) {
       nd <- ntile <- ensure2vector(ndminX)
       if(verbose)
         cat(paste("nd and ntile default to", nd[1], "*", nd[2], "\n"))
     } else {
       # find suitable divisors of the number of pixels
       nd <- ntile <- min2div(Nmin, ndminX, Nmax)
       if(any(nd >= Nmin)) { # none suitable
         if(verbose)
           cat("No suitable divisor of pixel dimensions\n")
         nd <- ntile <- ndminX
       }
     }
  } else {
    # both nd, ntile were given
    if(any(ntile > nd))
      warning("the number of tiles (ntile) exceeds the number of dummy points (nd)")
  }

  if(!ntile.given && quasi) {
    if(verbose) cat("Adjusting ntile because quasi=TRUE\n")
    ntile <- maxdiv(ntile, if(pixels) 2 else 1)
  } 
 
  if(!npix.given && pixels) 
    npix <- min2mul(nd, Nmin, Nmax)

  if(verbose) {
    if(!quasi)
      cat(paste("dummy points:",
                paste0(if(random) "stratified random in" else NULL,
                       "grid"),
                nd[1], "x", nd[2], "\n"))
    else
      cat(paste("dummy points:",
                nd[1], "x", nd[2], "=", prod(nd),
                "quasirandom points\n"))
    cat(paste("weighting tiles", ntile[1], "x", ntile[2], "\n"))
    if(pixels) cat(paste("pixel grid", npix[1], "x", npix[2], "\n"))
  }

  if(pixels) 
    return(list(nd=nd, ntile=ntile, npix=npix))
  else
    return(list(nd=nd, ntile=ntile, npix=npix))
}

  default.n.tiling
})


