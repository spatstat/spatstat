#'
#'     isclose.R
#'
#'    Determine whether each point has a close neighbour
#'
#'    $Revision: 1.7 $  $Date: 2016/10/15 09:25:25 $

is.close <- function(X, r, Y=NULL, ...) {
  UseMethod("is.close")
}

is.close.default <- function(X, r, Y=NULL, ..., periodic=FALSE) {
  trap.extra.arguments(...)
  if(!periodic) {
    nd <- if(is.null(Y)) nndist(X) else nncross(X, Y, what="dist")
    return(nd <= r)
  }
  if(is.null(Y)) {
    pd <- pairdist(X, periodic=TRUE)
    diag(pd) <- Inf
  } else {
    pd <- crossdist(X, Y, periodic=TRUE)
  }
  return(apply(pd <= r, 1, any))
}

is.close.ppp <- function(X, r, Y=NULL, ..., periodic=FALSE, sorted=FALSE) {
  trap.extra.arguments(...)
  nX <- npoints(X)
  if(nX <= 1) return(logical(nX))
  #' sort by increasing x coordinate
  cX <- coords(X)
  if(!sorted) {
    oo <- order(cX$x)
    cX <- cX[oo,,drop=FALSE]
  }
  if(is.null(Y)) {
    if(!periodic) {
      zz <- .C("isXclose",
               n = as.integer(nX),
               x = as.double(cX$x),
               y = as.double(cX$y),
               r = as.double(r),
               t = as.integer(integer(nX)))
    } else {
      b <- sidelengths(Frame(X))
      zz <- .C("isXpclose",
               n = as.integer(nX),
               x = as.double(cX$x),
               y = as.double(cX$y),
               r = as.double(r),
               b = as.double(b),
               t = as.integer(integer(nX)))
    }
  } else {
    stopifnot(is.ppp(Y))
    nY <- npoints(Y)
    if(nY == 0) return(logical(nX))
    cY <- coords(Y)
    #' sort Y by increasing x coordinate
    if(!sorted) {
      ooY <- order(cY$x)
      cY <- cY[ooY, , drop=FALSE]
    }
    if(!periodic) {
      zz <- .C("isXYclose",
               n1 = as.integer(nX),
               x1 = as.double(cX$x),
               y1 = as.double(cX$y),
               n2 = as.integer(nY),
               x2 = as.double(cY$x),
               y2 = as.double(cY$y),
               r = as.double(r),
               t = as.integer(integer(nX)))
    } else {
      bX <- sidelengths(Frame(X))
      bY <- sidelengths(Frame(Y))
      if(any(bX != bY))
        warning("Windows are not equal: periodic distance may be erroneous")
      zz <- .C("isXYpclose",
               n1 = as.integer(nX),
               x1 = as.double(cX$x),
               y1 = as.double(cX$y),
               n2 = as.integer(nY),
               x2 = as.double(cY$x),
               y2 = as.double(cY$y),
               r = as.double(r),
               b = as.double(bX),
               t = as.integer(integer(nX)))
    }
  }
  tt <- as.logical(zz$t)
  if(sorted) return(tt)
  #' reinstate original order
  ans <- logical(nX)
  ans[oo] <- tt
  return(ans)
}

is.close.pp3 <- function(X, r, Y=NULL, ..., periodic=FALSE, sorted=FALSE) {
  trap.extra.arguments(...)
  nX <- npoints(X)
  if(nX <= 1) return(logical(nX))
  cX <- coords(X)
  if(!sorted) {
    #' sort by increasing x coordinate
    oo <- order(cX$x)
    cX <- cX[oo,,drop=FALSE]
  }
  if(is.null(Y)) {
    if(!periodic) {
      zz <- .C("isX3close",
               n = as.integer(nX),
               x = as.double(cX$x),
               y = as.double(cX$y),
               z = as.double(cX$z),
               r = as.double(r),
               t = as.integer(integer(nX)))
    } else {
      b <- sidelengths(as.box3(X))
      zz <- .C("isX3pclose",
               n = as.integer(nX),
               x = as.double(cX$x),
               y = as.double(cX$y),
               z = as.double(cX$z),
               r = as.double(r),
               b = as.double(b), 
               t = as.integer(integer(nX)))
    }
  } else {
    stopifnot(is.pp3(Y))
    nY <- npoints(Y)
    if(nY == 0) return(logical(nX))
    cY <- coords(Y)
    if(!sorted) {
      #' sort Y by increasing x coordinate
      ooY <- order(cY$x)
      cY <- cY[ooY, , drop=FALSE]
    }
    if(!periodic) {
      zz <- .C("isXYclose",
               n1 = as.integer(nX),
               x1 = as.double(cX$x),
               y1 = as.double(cX$y),
               z1 = as.double(cX$z),
               n2 = as.integer(nY),
               x2 = as.double(cY$x),
               y2 = as.double(cY$y),
               z2 = as.double(cY$z),
               r = as.double(r),
               t = as.integer(integer(nX)))
    } else {
      bX <- sidelengths(as.box3(X))
      bY <- sidelengths(as.box3(Y))
      if(any(bX != bY))
        warning("Domains are not equal: periodic distance may be erroneous")
      zz <- .C("isXYclose",
               n1 = as.integer(nX),
               x1 = as.double(cX$x),
               y1 = as.double(cX$y),
               z1 = as.double(cX$z),
               n2 = as.integer(nY),
               x2 = as.double(cY$x),
               y2 = as.double(cY$y),
               z2 = as.double(cY$z),
               r = as.double(r),
               b = as.double(bX),
               t = as.integer(integer(nX)))
    }
  }
  tt <- as.logical(zz$t)
  if(sorted) return(tt)
  #' reinstate original order
  ans <- logical(nX)
  ans[oo] <- tt
  return(ans)
}


  
