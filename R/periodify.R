#
# periodify.R
#
# replicate a pattern periodically
#
#  $Revision: 1.3 $  $Date: 2011/04/17 05:52:50 $
#

periodify <- function(X, ...) {
  UseMethod("periodify")
}

periodify.ppp <- function(X, nx=1, ny=1, ...,
                          combine=TRUE, warn=TRUE, check=TRUE, 
                          ix=(-nx):nx, iy=(-ny):ny,
                          ixy=expand.grid(ix=ix,iy=iy)) {
  # sanity checks
  if(!missing(nx) || !missing(ny)) {
    if(is.null(nx)) nx <- 1
    if(is.null(ny)) ny <- 1
    if(length(nx) != 1 || length(ny) != 1)
      stop("nx and ny should be single integers")
    if(nx != round(nx) || ny != round(ny))
      stop("nx and ny should be integers")
  }
  force(ixy)
  W <- X$window
  isrect <- (W$type == "rectangle")
  if(warn && combine && !isrect)
    warning("X has a non-rectangular window")
  else 
   isrect <- isrect && all(diff(nx) == 1) && all(diff(ny) == 1)
  width <- diff(W$xrange)
  height <- diff(W$yrange)
  shifts <- cbind(ixy[,1] * width, ixy[,2] * height)
  Xshift <- list()
  for(i in 1:nrow(shifts))
    Xshift[[i]] <- shift(X, vec=as.numeric(shifts[i, ]))
  if(!combine)
    return(Xshift)
  Wnew <- if(isrect) {
    owin(range(range(W$xrange) + range(shifts[,1])),
         range(range(W$yrange) + range(shifts[,2])))
  } else NULL
  Z <- do.call(superimpose, append(Xshift, list(W=Wnew, check=check)))
  return(Z)
}

periodify.psp <- function(X, nx=1, ny=1, ...,
                          combine=TRUE, warn=TRUE, check=TRUE,
                          ix=(-nx):nx, iy=(-ny):ny,
                          ixy=expand.grid(ix=ix,iy=iy)) {
  # sanity checks
  if(!missing(nx) || !missing(ny)) {
    if(is.null(nx)) nx <- 1
    if(is.null(ny)) ny <- 1
    if(length(nx) != 1 || length(ny) != 1)
      stop("nx and ny should be single integers")
    if(nx != round(nx) || ny != round(ny))
      stop("nx and ny should be integers")
  }
  force(ixy)
  W <- X$window
  isrect <- (W$type == "rectangle")
  if(warn && combine && !isrect)
    warning("X has a non-rectangular window")
  else 
   isrect <- isrect && all(diff(nx) == 1) && all(diff(ny) == 1)
  width <- diff(W$xrange)
  height <- diff(W$yrange)
  shifts <- cbind(ixy[,1] * width, ixy[,2] * height)
  Xshift <- list()
  for(i in 1:nrow(shifts))
    Xshift[[i]] <- shift(X, vec=as.numeric(shifts[i, ]))
  if(!combine)
    return(Xshift)
  Wnew <- if(isrect) {
    owin(range(range(W$xrange) + range(shifts[,1])),
         range(range(W$yrange) + range(shifts[,2])))
  } else NULL
  Z <- do.call(superimpose, append(Xshift, list(W=Wnew, check=check)))
  return(Z)
}

periodify.owin <- function(X, nx=1, ny=1, ...,
                          combine=TRUE, warn=TRUE,
                          ix=(-nx):nx, iy=(-ny):ny,
                          ixy=expand.grid(ix=ix,iy=iy)) {
  # sanity checks
  if(!missing(nx) || !missing(ny)) {
    if(is.null(nx)) nx <- 1
    if(is.null(ny)) ny <- 1
    if(length(nx) != 1 || length(ny) != 1)
      stop("nx and ny should be single integers")
    if(nx != round(nx) || ny != round(ny))
      stop("nx and ny should be integers")
  }
  force(ixy)
  isrect <- (X$type == "rectangle")
  if(warn && combine && !isrect)
    warning("X is not rectangular")
  else 
    isrect <- isrect && all(diff(nx) == 1) && all(diff(ny) == 1)
  width <- diff(X$xrange)
  height <- diff(X$yrange)
  shifts <- cbind(ixy[,1] * width, ixy[,2] * height)
  if(combine) {
    if(isrect) {
      # result is a rectangle
      Y <-  owin(range(range(X$xrange) + range(shifts[,1])),
                    range(range(X$yrange) + range(shifts[,2])))
    } else {
      # result is another type of window
      for(i in 1:nrow(shifts)) {
        Xi <- shift(X, vec=as.numeric(shifts[i, ]))
        Y <- if(i == 1) Xi else union.owin(Y, Xi)
      }
    }
  } else {
    # result is a list
    Y <- list()
    for(i in 1:nrow(shifts))
      Y[[i]] <- shift(X, vec=as.numeric(shifts[i, ]))
  }
  return(Y)
}

