##
##    parameters.R
##
##   $Revision: 1.2 $ $Date: 2015/05/08 04:27:15 $
##

parameters <- function(model, ...) {
  UseMethod("parameters")
}

parameters.ppm <- function(model, ...) {
  ss <- summary(model, quick="no variances")
  out <- c(list(trend=ss$trend$value),
           ss$covfunargs,
           ss$interaction$interaction$par,
           ss$interaction$sensible$param)
  return(out)
}

parameters.kppm <- function(model, ...) {
  ss <- summary(model, quick="no variances")
  out <- c(list(trend=ss$trend$trend$value),
           ss$covfunargs,
           ss$clustpar,
           ss$clustargs,
           list(mu=ss$mu))
  return(out)
}


