#
#
#   rescale.R
#
#   $Revision: 1.8 $ $Date: 2019/02/10 06:42:26 $
#
#

rescale <- function(X, s, unitname) {
  UseMethod("rescale")
}

rescale.ppp <- function(X, s, unitname) {
  if(missing(unitname)) unitname <- NULL
  if(missing(s) || is.null(s)) s <- 1/unitname(X)$multiplier
  Y <- affine.ppp(X, mat=diag(c(1/s,1/s)), rescue=FALSE)
  unitname(Y) <- rescale(unitname(X), s, unitname)
  return(Y)
}

rescale.owin <- function(X, s, unitname) {
  if(missing(unitname)) unitname <- NULL
  if(missing(s) || is.null(s)) s <- 1/unitname(X)$multiplier
  Y <- affine.owin(X, mat=diag(c(1/s,1/s)), rescue=FALSE)
  unitname(Y) <- rescale(unitname(X), s, unitname)
  return(Y)
}

rescale.im <- function(X, s, unitname) {
  if(missing(unitname)) unitname <- NULL
  if(missing(s) || is.null(s)) s <- 1/unitname(X)$multiplier
  Y <- X
  Y$xrange <- X$xrange/s
  Y$yrange <- X$yrange/s
  Y$xstep  <- X$xstep/s
  Y$ystep  <- X$ystep/s
  Y$xcol   <- X$xcol/s
  Y$yrow   <- X$yrow/s
  unitname(Y) <- rescale(unitname(X), s, unitname)
  return(Y)
}

rescale.psp <- function(X, s, unitname) {
  if(missing(unitname)) unitname <- NULL
  if(missing(s) || is.null(s)) s <- 1/unitname(X)$multiplier
  Y <- affine.psp(X, mat=diag(c(1/s,1/s)), rescue=FALSE)
  unitname(Y) <- rescale(unitname(X), s, unitname)
  return(Y)
}
  
rescale.unitname <- function(X, s, unitname) {
  if(!missing(unitname) && !is.null(unitname)) return(as.unitname(unitname))
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

