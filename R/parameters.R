##
##    parameters.R
##
##   $Revision: 1.1 $ $Date: 2015/04/25 21:52:59 $
##

parameters <- function(model, ...) {
  UseMethod("parameters")
}

parameters.ppm <- function(model, ...) {
  ss <- summary(model)
  out <- c(list(trend=ss$trend$value),
           ss$covfunargs,
           ss$interaction$interaction$par,
           ss$interaction$sensible$param)
  return(out)
}

parameters.kppm <- function(model, ...) {
  ss <- summary(model)
  out <- c(list(trend=ss$trend$trend$value),
           ss$covfunargs,
           ss$clustpar,
           ss$clustargs,
           list(mu=ss$mu))
  return(out)
}


