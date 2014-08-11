#
#
#   rescale.R
#
#   $Revision: 1.5 $ $Date: 2013/04/25 06:37:43 $
#
#

rescale <- function(X, s) {
  UseMethod("rescale")
}

rescale.ppp <- function(X, s) {
  if(missing(s)) s <- 1/unitname(X)$multiplier
  Y <- affine.ppp(X, mat=diag(c(1/s,1/s)))
  unitname(Y) <- rescale(unitname(X), s)
  return(Y)
}

rescale.owin <- function(X, s) {
  if(missing(s)) s <- 1/unitname(X)$multiplier
  Y <- affine.owin(X, mat=diag(c(1/s,1/s)))
  unitname(Y) <- rescale(unitname(X), s)
  return(Y)
}

rescale.im <- function(X, s) {
  if(missing(s)) s <- 1/unitname(X)$multiplier
  Y <- X
  Y$xrange <- X$xrange/s
  Y$yrange <- X$yrange/s
  Y$xstep  <- X$xstep/s
  Y$ystep  <- X$ystep/s
  Y$xcol   <- X$xcol/s
  Y$yrow   <- X$yrow/s
  unitname(Y) <- rescale(unitname(X), s)
  return(Y)
}

rescale.psp <- function(X, s) {
  if(missing(s)) s <- 1/unitname(X)$multiplier
  Y <- affine.psp(X, mat=diag(c(1/s,1/s)))
  unitname(Y) <- rescale(unitname(X), s)
  return(Y)
}
  
rescale.units <- function(X, s) {
  if(summary(X)$vanilla)
    return(X)
  if(missing(s)) {
    X$multiplier <- 1
  } else {
    if(!is.numeric(s) || length(s) != 1 || s <= 0)
      stop("s should be a positive number")
    X$multiplier <- s * X$multiplier
  }
  return(X)
}


