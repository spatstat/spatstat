#
#	distbdry.S		Distance to boundary
#
#	$Revision: 4.36 $	$Date: 2013/04/25 06:37:43 $
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
                 result <- rep.int(Inf, X$n)
                 bdry <- window$bdry
                 for(i in seq_along(bdry)) {
                   polly <- bdry[[i]]
                   px <- polly$x
                   py <- polly$y
                   nsegs <- length(px)
                   for(j in seq_len(nsegs)) {
                     j1 <- if(j < nsegs) j + 1 else 1
                     seg <- c(px[j],  py[j],
                              px[j1], py[j1])
                     result <- pmin.int(result, distppl(xy, seg))
                   }
                 }
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

"bdist.pixels" <- function(w, ..., style="image") {
	verifyclass(w, "owin")

        masque <- as.mask(w, ...)
        
        switch(w$type,
               mask = {
                 neg <- complement.owin(masque)
                 m <- exactPdt(neg)
                 b <- pmin.int(m$d,m$b)
               },
               rectangle = {
                 x <- raster.x(masque)
                 y <- raster.y(masque)
                 xmin <- w$xrange[1]
                 xmax <- w$xrange[2]
                 ymin <- w$yrange[1]
                 ymax <- w$yrange[2]
                 b <- pmin.int(x - xmin, xmax - x, y - ymin, ymax - y)
               },
               polygonal = {
                 # set up pixel raster
                 x <- as.vector(raster.x(masque))
                 y <- as.vector(raster.y(masque))
                 b <- numeric(length(x))
                 # test each pixel in/out, analytically
                 inside <- inside.owin(x, y, w)
                 # compute distances for these pixels
                 xy <- cbind(x[inside], y[inside])
                 dxy <- rep.int(Inf, sum(inside))
                 bdry <- w$bdry
                 for(i in seq_along(bdry)) {
                   polly <- bdry[[i]]
                   nsegs <- length(polly$x)
                   for(j in 1:nsegs) {
                     j1 <- if(j < nsegs) j + 1 else 1
                     seg <- c(polly$x[j],  polly$y[j],
                              polly$x[j1], polly$y[j1])
                     dxy <- pmin.int(dxy, distppl(xy, seg))
                   }
                 }
                 b[inside] <- dxy
               },
               stop("unrecognised window type", w$type)
               )

        # reshape it
        b <- matrix(b, nrow=masque$dim[1], ncol=masque$dim[2])

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
  if(!is.numeric(r) || length(r) != 1)
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
  switch(w$type,
         rectangle={
           return(owin(xr, yr,
                       poly=list(x=w$xrange[c(1,2,2,1)],
                                 y=w$yrange[c(1,1,2,2)])))
         },
         polygonal={
           return(owin(xr, yr, poly=w$bdry, check=FALSE))
         },
         mask={
           newseq <- function(oldseq, newrange, dstep) {
             oldends <- range(oldseq)
             nleft <- max(0, floor((oldends[1] - newrange[1])/dstep))
             nright <- max(0, floor((newrange[2] - oldends[2])/dstep))
             newstart <- max(oldends[1] - nleft * dstep, newrange[1])
             newend <- min(oldends[2] + nright * dstep, newrange[2])
             seq(from=newstart, by=dstep, to=newend)
           }
           xcol <- newseq(w$xcol, xr, mean(diff(w$xcol)))
           yrow <- newseq(w$yrow, yr, mean(diff(w$yrow)))
           newmask <- as.mask(xy=list(x=xcol, y=yrow))
           xx <- raster.x(newmask)
           yy <- raster.y(newmask)
           newmask$m <- inside.owin(xx, yy, w)
           return(newmask)
         }
         )
}
  
    
  
