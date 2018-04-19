#'
#'   dffit.R
#'
#'   $Revision: 1.1 $ $Date: 2018/04/19 05:04:59 $
#'
#'   Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2018

dffit <- function(object, ...) UseMethod("dffit")

dffit.ppm <- function(object, ..., collapse=FALSE, dfb=NULL) {
  db <- dfb %orifnull% dfbetas(object, ...)
  Z <- model.matrix(object, drop=FALSE, irregular=TRUE)
  if(!all(dim(db) == dim(Z))) {
    #' internal error - mismatch in quadrature points - fall back
    U <- db$loc
    Z <- sapply(model.images(object, irregular=TRUE), "[", i=U)
  }
  #' ensure 0 * (-Inf) = 0
  if(any(a <- (db$val == 0) & (Z == -Inf))) 
    Z[a] <- 0
  #' smoothed density must be handled separately
  sm <- attr(db, "smoothdensity")
  attr(db, "smoothdensity") <- NULL
  #' do the main calculation
  Y <- db * Z
  #' also calculate the smoothed density
  if(!is.null(sm)) {
    ZZ <- model.images(object, irregular=TRUE)
    HH <- mapply(harmonise, ZZ=ZZ, sm=sm)
    sm <- mapply("*", e1=HH$sm, e2=HH$ZZ, SIMPLIFY=FALSE)
    attr(Y, "smoothdensity") <- as.solist(sm)
  }
  if(collapse) Y <- Reduce("+", unstack(Y))
  return(Y)
}

  
