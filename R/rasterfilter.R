#'
#'    rasterfilter.R
#'
#'    raster filters implemented directly
#'
#'    $Revision: 1.5 $ $Date: 2017/11/18 07:17:18 $
#'

rasterfilter <- function(X, f) {
  X <- as.im(X)
  dimX <- dim(X)
  f <- as.matrix(f)
  if(!all(dim(f) == 3))
    stop("f should be a 3 x 3 matrix or image")
  #' handle NA
  v <- as.double(X$v)
  if(hasna <- anyNA(v)) {
    isna <- is.na(v)
    v[isna] <- 0
  }
  #' compute
  z <- .C("raster3filter",
          nx = as.integer(dimX[2]),
          ny = as.integer(dimX[1]),
          a  = as.double(v),
          w  = as.double(f),
          b  = as.double(numeric(prod(dimX))),
          PACKAGE="spatstat")
  z <- z$b
  #' handle NA
  if(hasna)
    z[isna] <- NA
  # replace
  X[] <- z
  return(X)
}

#'  antialiasing
smudge <- function(X) {
  stopifnot(is.im(X))
  xstep <- X$xstep
  ystep <- X$ystep
  #' choose a very small bandwidth
  sigma <- min(xstep, ystep)/2
  #' match variance: 2 p step^2 = sigma^2
  px <- sigma^2/(2 * xstep^2)
  py <- sigma^2/(2 * ystep^2)
  f <- outer(c(py, 1-2*py, py), c(px, 1-2*px, px), "*")
  #' compute 
  Z <- rasterfilter(X, f)
  attr(Z, "sigma") <- sigma
  return(Z)
}
