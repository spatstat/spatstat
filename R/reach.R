#
#   reach.R
#
#  $Revision: 1.8 $   $Date: 2007/10/24 09:41:15 $
#

reach <- function(x, ...) {
  UseMethod("reach")
}

reach.interact <- function(x, ...) {
  verifyclass(x, "interact")
  irange <- x$irange
  if(is.null(irange))
    return(Inf)
  if(!is.function(irange))
    stop("Internal error - x$irange is not a function")
  ir <- irange(x)
  if(is.na(ir))
    ir <- Inf
  return(ir)
}

reach.ppm <- function(x, ..., epsilon=0) {
  verifyclass(x, "ppm")
  
  # Poisson case
  if(is.poisson.ppm(x))
    return(0)

  # extract info
  inte <- x$interaction
  coeffs <- coef(x)

  if(newstyle.coeff.handling(inte)) {
    # extract only interaction coefficients
    Vnames <- x$internal$Vnames
    coeffs <- coeffs[Vnames]
  } 
  
  # apply 'irange' function
  irange <- inte$irange
  if(is.null(irange))
    return(Inf)
  ir <- irange(inte, coeffs, epsilon=epsilon)

  if(is.na(ir))
    ir <- Inf

  return(ir)
}



