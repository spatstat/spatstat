#
#	distbdry.S		Distance to boundary
#
#	$Revision: 4.45 $	$Date: 2019/02/06 10:53:48 $
#
# -------- functions ----------------------------------------
#
#	bdist.points()
#                       compute vector of distances 
#			from each point of point pattern
#                       to boundary of window
#
#       bdist.pixels()
#                       compute matrix of distances from each pixel
#                       to boundary of window
#
#       erodemask()    erode the window mask by a distance r
#                       [yields a new window]
#
#
# 
"bdist.points"<-
function(X)
{
	verifyclass(X, "ppp") 
        if(X$n == 0)
          return(numeric(0))
	x <- X$x
	y <- X$y
	window <- X$window
        switch(window$type,
               rectangle = {
		xmin <- min(window$xrange)
		xmax <- max(window$xrange)
		ymin <- min(window$yrange)
		ymax <- max(window$yrange)
		result <- pmin.int(x - xmin, xmax - x, y - ymin, ymax - y)
               },
               polygonal = {
                 xy <- cbind(x,y)
                 ll <- edges(window)$ends
                 result <- distppllmin(xy, ll)$min.d
               },
               mask = {
                 b <- bdist.pixels(window, style="matrix")
                 loc <- nearest.raster.point(x,y,window)
                 result <- b[cbind(loc$row, loc$col)]
               },
               stop("Unrecognised window type", window$type)
               )
        return(result)
}

"bdist.pixels" <- function(w, ..., style="image",
                           method=c("C", "interpreted")) {
	verifyclass(w, "owin")

        masque <- as.mask(w, ...)
        
        switch(w$type,
               mask = {
                 neg <- complement.owin(masque)
                 m <- exactPdt(neg)
                 b <- pmin.int(m$d,m$b)
               },
               rectangle = {
                 rxy <- rasterxy.mask(masque)
                 x <- rxy$x
                 y <- rxy$y
                 xmin <- w$xrange[1L]
                 xmax <- w$xrange[2L]
                 ymin <- w$yrange[1L]
                 ymax <- w$yrange[2L]
                 b <- pmin.int(x - xmin, xmax - x, y - ymin, ymax - y)
               },
               polygonal = {
                 # set up pixel raster
                 method <- match.arg(method)
                 rxy <- rasterxy.mask(masque)
                 x <- rxy$x
                 y <- rxy$y
                 b <- numeric(length(x))
                 # test each pixel in/out, analytically
                 inside <- inside.owin(x, y, w)
                 # compute distances for these pixels
                 xy <- cbind(x[inside], y[inside])
                 switch(method,
                        C = {
                          #' C code
                          ll <- as.data.frame(edges(w))
                          dxy <- distppllmin(xy, ll)$min.d
                        },
                        interpreted = {
                          #' ancient R code
                          dxy <- rep.int(Inf, sum(inside))
                          bdry <- w$bdry
                          for(i in seq_along(bdry)) {
                            polly <- bdry[[i]]
                            nsegs <- length(polly$x)
                            for(j in 1:nsegs) {
                              j1 <- if(j < nsegs) j + 1L else 1L
                              seg <- c(polly$x[j],  polly$y[j],
                                       polly$x[j1], polly$y[j1])
                              dxy <- pmin.int(dxy, distppl(xy, seg))
                            }
                          }
                        })
                 b[inside] <- dxy
               },
               stop("unrecognised window type", w$type)
               )

        # reshape it
        b <- matrix(b, nrow=masque$dim[1L], ncol=masque$dim[2L])

        switch(style,
               coords={
                 # format which can be plotted by image(), persp() etc
                 return(list(x=masque$xcol, y=masque$yrow, z=t(b)))
               },
               matrix={
                 # return matrix (for internal use by package)
                 return(b)
               },
               image={
                 bim <- im(b, xcol=masque$xcol, yrow=masque$yrow,
                           unitname=unitname(masque))
                 return(bim)
               },
               stop(paste("Unrecognised option for style:", style)))
} 

erodemask <- function(w, r, strict=FALSE) {
  # erode a binary image mask without changing any other entries
  verifyclass(w, "owin")
  if(w$type != "mask")
    stop(paste("window w is not of type", sQuote("mask")))
  if(!is.numeric(r) || length(r) != 1L)
    stop("r must be a single number")
  if(r < 0)
    stop("r must be nonnegative")
        
  bb <- bdist.pixels(w, style="matrix")

  if(r > max(bb))
    warning("eroded mask is empty")

  if(identical(strict, TRUE))
    w$m <- (bb > r)
  else 
    w$m <- (bb >= r)
  return(w)
}

"Frame<-.owin" <- function(X, value) {
  stopifnot(is.rectangle(value))
  W <- Frame(X)
  if(!is.subset.owin(W, value))
    W <- intersect.owin(W, value)
  rebound.owin(X, value)
}

rebound.owin <- local({

  rebound.owin <- function(x, rect) {
    w <- x
    verifyclass(rect, "owin")
    if(is.empty(w))
      return(emptywindow(rect))
    verifyclass(w, "owin")
    if(!is.subset.owin(as.rectangle(w), rect)) {
      bb <- boundingbox(w)
      if(!is.subset.owin(bb, rect))
        stop(paste("The new rectangle",
                   sQuote("rect"),
                   "does not contain the window",
                   sQuote("win")))
    }
    xr <- rect$xrange
    yr <- rect$yrange
    ## determine unitname
    uu <- list(unitname(x), unitname(rect))
    uu <- unique(uu[sapply(uu, is.vanilla)])
    if(length(uu) > 1) {
      warning("Arguments of rebound.owin have incompatible unitnames",
              call.=FALSE)
      uu <- list()
    }
    un <- if(length(uu)) uu[[1]] else NULL
    ## 
    switch(w$type,
           rectangle={
             return(owin(xr, yr,
                         poly=list(x=w$xrange[c(1L,2L,2L,1L)],
                                   y=w$yrange[c(1L,1L,2L,2L)]),
                         unitname = un,
                         check=FALSE))
           },
           polygonal={
             return(owin(xr, yr, poly=w$bdry, unitname=un, check=FALSE))
           },
           mask={
             xcol <- newseq(w$xcol, xr)
             yrow <- newseq(w$yrow, yr)
             newmask <- as.mask(xy=list(x=xcol, y=yrow))
             xx <- rasterx.mask(newmask)
             yy <- rastery.mask(newmask)
             newmask$m <- inside.owin(xx, yy, w)
             unitname(newmask) <- un
             return(newmask)
           }
           )
  }

  newseq <- function(oldseq, newrange) {
    oldrange <- range(oldseq)
    dstep <- mean(diff(oldseq))
    nleft <- max(0, floor((oldrange[1L] - newrange[1L])/dstep))
    nright <- max(0, floor((newrange[2L] - oldrange[2L])/dstep))
    newstart <- max(oldrange[1L] - nleft * dstep, newrange[1L])
    newend <- min(oldrange[2L] + nright * dstep, newrange[2L])
    seq(from=newstart, by=dstep, to=newend)
  }

  rebound.owin
})

    
  
