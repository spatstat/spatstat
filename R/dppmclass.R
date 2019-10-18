#' Methods and support for class dppm
#'
#' $Revision: 1.5 $ $Date: 2019/10/18 03:41:43 $

is.dppm <- function(x) { inherits(x, "dppm") }

plot.dppm <- function (x, ..., what = c("intensity", "statistic")){
    objectname <- short.deparse(substitute(x))
    if(missing(what) && is.stationary(x))
        what <- "statistic"
    plot.kppm(x, ..., xname = objectname, what = what)
}

Kmodel.dppm <- function (model, ...){
    Kmodel(model$fitted, W=model$window)
}

pcfmodel.dppm <- function (model, ...){
    pcfmodel(model$fitted, W=model$window)
}

intensity.dppm <- function (X, ...){
    return(intensity(X$fitted))
}

reach.dppm <- function(x, ...){
    reach(x$fitted, ...)
}

repul <- function(model, ...) {
  UseMethod("repul")
}

repul.dppm <- function(model, ...) {
  g <- pcfmodel(model)
  f <- function(x) { 2 * pi * x * (1 - g(x)) }
  rmax <- reach(model)
  h <- integrate(f, 0, rmax)$value
  lam <- intensity(model)
  ans <- h * lam
  return(ans)
}

#' print.dppm is identical to print.kppm and defined in kppm.R
#' summary.dppm is defined in summary.dppm.R
