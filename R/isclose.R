#'
#'     isclose.R
#'
#'    Determine whether each point has a close neighbour
#'
#'    $Revision: 1.4 $  $Date: 2016/10/14 07:36:22 $

is.close <- function(X, r, Y=NULL) {
  UseMethod("is.close")
}

is.close.default <- function(X, r, Y=NULL) {
  nd <- if(is.null(Y)) nndist(X) else nncross(X, Y, what="dist")
  return(nd <= r)
}

is.close.ppp <- function(X, r, Y=NULL) {
  nX <- npoints(X)
  if(nX <= 1) return(logical(nX))
  #' sort by increasing x coordinate
  cX <- coords(X)
  oo <- order(cX$x)
  cX <- cX[oo,,drop=FALSE]
  if(is.null(Y)) {
    zz <- .C("isXclose",
             n = as.integer(nX),
             x = as.double(cX$x),
             y = as.double(cX$y),
             r = as.double(r),
             t = as.integer(integer(nX)))
    tt <- as.logical(zz$t)
  } else {
    stopifnot(is.ppp(Y))
    #' sort Y by increasing x coordinate
    cY <- coords(Y)
    ooY <- order(cY$x)
    cY <- cY[ooY, , drop=FALSE]
    nY <- npoints(Y)
    zz <- .C("isXYclose",
             n1 = as.integer(nX),
             x1 = as.double(cX$x),
             y1 = as.double(cX$y),
             n2 = as.integer(nY),
             x2 = as.double(cY$x),
             y2 = as.double(cY$y),
             r = as.double(r),
             t = as.integer(integer(nX)))
    tt <- as.logical(zz$t)
  }
  #' reinstate original order
  ans <- logical(nX)
  ans[oo] <- tt
  return(ans)
}

is.close.pp3 <- function(X, r, Y=NULL) {
  nX <- npoints(X)
  if(nX <= 1) return(logical(nX))
  #' sort by increasing x coordinate
  cX <- coords(X)
  oo <- order(cX$x)
  cX <- cX[oo,,drop=FALSE]
  if(is.null(Y)) {
    zz <- .C("isX3close",
             n = as.integer(nX),
             x = as.double(cX$x),
             y = as.double(cX$y),
             z = as.double(cX$z),
             r = as.double(r),
             t = as.integer(integer(nX)))
    tt <- as.logical(zz$t)
  } else {
    stopifnot(is.pp3(Y))
    #' sort Y by increasing x coordinate
    cY <- coords(Y)
    ooY <- order(cY$x)
    cY <- cY[ooY, , drop=FALSE]
    nY <- npoints(Y)
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
    tt <- as.logical(zz$t)
  }
  #' reinstate original order
  ans <- logical(nX)
  ans[oo] <- tt
  return(ans)
}


  
