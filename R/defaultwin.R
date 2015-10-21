#
#
#  defaultwin.R
#
#   $Revision: 1.10 $   $Date: 2015/10/21 09:06:57 $
#

default.expand <- function(object, m=2, epsilon=1e-6, w=Window(object)) {
  stopifnot(is.ppm(object) || inherits(object, "rmhmodel"))
  # no expansion necessary if model is Poisson
  if(is.poisson(object))
    return(.no.expansion)
  # default is no expansion if model is nonstationary
  if(!is.stationary(object))
    return(.no.expansion)
  
# Redundant since a non-expandable model is non-stationary
#  if(!is.expandable(object))
#    return(.no.expansion)
  
  # rule is to expand data window by distance d = m * reach
  rr <- reach(object, epsilon=epsilon)
  if(!is.finite(rr))
    return(rmhexpand())
  if(!is.numeric(m) || length(m) != 1 || m < 1)
    stop("m should be a single number >= 1")
  mr <- m * rr
  rule <- rmhexpand(distance = mr)
  # 
  if(is.owin(w)) {
    # apply rule to window
    wplus <- expand.owin(w, rule)
    # save as new expansion rule
    rule <- rmhexpand(wplus)
  }
  return(rule)
}

default.clipwindow <- function(object, epsilon=1e-6) {
  stopifnot(is.ppm(object) || inherits(object, "rmhmodel"))
  # data window
  w <- as.owin(object)
  if(is.null(w)) return(NULL)
  # interaction range of model
  rr <- reach(object, epsilon=epsilon)
  if(!is.finite(rr))
    return(NULL)
  if(rr == 0)
    return(w)
  else
    return(erosion(w, rr))
}

  
