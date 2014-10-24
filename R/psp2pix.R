#
# psp2pix.R
#
#  $Revision: 1.7 $  $Date: 2014/10/24 00:22:30 $
#
#

as.mask.psp <- function(x, W=NULL, ...) {
  L <- as.psp(x)
  if(is.null(W))
    W <- as.owin(L)
  else
    W <- as.owin(W)
  W <- as.mask(W, ...)

  ends <- L$ends
  nseg <- nrow(ends)
  
  if(nseg == 0) {
    # empty
    W$m[] <- FALSE
    return(W)
  }
    
  x0 <- (ends$x0 - W$xrange[1])/W$xstep
  x1 <- (ends$x1 - W$xrange[1])/W$xstep
  y0 <- (ends$y0 - W$yrange[1])/W$ystep
  y1 <- (ends$y1 - W$yrange[1])/W$ystep
  nr <- W$dim[1]
  nc <- W$dim[2]
  zz <- .C("seg2pixI",
           ns=as.integer(nseg),
           x0=as.double(x0),
           y0=as.double(y0),
           x1=as.double(x1),
           y1=as.double(y1),
           nx=as.integer(nc),
           ny=as.integer(nr),
           out=as.integer(integer(nr * nc)))
  mm <- matrix(zz$out, nr, nc)
  # intersect with existing window
  W$m <- W$m & mm
  W
}


pixellate.psp <- function(x, W=NULL, ..., weights=NULL) {
  L <- as.psp(x)
  if(is.null(W))
    W <- as.owin(L)
  else
    W <- as.owin(W)
  W <- as.mask(W, ...)

  Z <- as.im(W)

  ends <- L$ends
  nseg <- nrow(ends)

  if(nseg == 0) {
    # empty
    Z$v[] <- 0
    return(Z)
  }
  
  if(is.null(weights))
    weights <- rep.int(1, nseg)
  else {
    if(!is.numeric(weights)) stop("weights must be numeric")
    if(any(is.na(weights))) stop("weights must not be NA")
    if(!all(is.finite(weights))) stop("weights must not be infinite")
    if(length(weights) == 1)
      weights <- rep.int(weights, nseg)
    else if(length(weights) != nseg)
      stop(paste("weights vector has length", length(weights),
                 "but there are", nseg, "line segments"))
  }
      
  x0 <- (ends$x0 - Z$xrange[1])/Z$xstep
  x1 <- (ends$x1 - Z$xrange[1])/Z$xstep
  y0 <- (ends$y0 - Z$yrange[1])/Z$ystep
  y1 <- (ends$y1 - Z$yrange[1])/Z$ystep
  nr <- Z$dim[1]
  nc <- Z$dim[2]
  zz <- .C("seg2pixL",
           ns=as.integer(nseg),
           x0=as.double(x0),
           y0=as.double(y0),
           x1=as.double(x1),
           y1=as.double(y1),
           weights=as.double(weights),
           pixwidth=as.double(Z$xstep),
           pixheight=as.double(Z$ystep),
           nx=as.integer(nc),
           ny=as.integer(nr),
           out=as.double(numeric(nr * nc)))
  mm <- matrix(zz$out, nr, nc)
  mm[is.na(Z$v)] <- NA
  # intersect with existing window
  Z$v <- mm
  Z
}


