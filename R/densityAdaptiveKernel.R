#'
#'   densityAdaptiveKernel.R
#'
#'   $Revision: 1.2 $  $Date: 2019/02/06 04:03:25 $
#'
#'   currently a stub!

densityAdaptiveKernel <- function(X, ...) {
  UseMethod("densityAdaptiveKernel")
}

densityAdaptiveKernel.ppp <- function(X, ...) {
  density.ppp(X, ...)
}
