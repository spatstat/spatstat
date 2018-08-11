#'
#'    hasenvelope.R
#'
#'    A simple class of objects which contain additional envelope data
#' 
#'    $Revision: 1.1 $ $Date: 2015/10/05 06:20:31 $

hasenvelope <- function(X, E=NULL) {
  if(inherits(E, "envelope")) {
    attr(X, "envelope") <- E
    class(X) <- c("hasenvelope", class(X))
  }
  return(X)
}

print.hasenvelope <- function(x, ...) {
  NextMethod("print")
  splat("[Object contains simulation envelope data]")
  return(invisible(NULL))
}

envelope.hasenvelope <- function(Y, ..., Yname=NULL) {
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  E <- attr(Y, "envelope")
  return(envelope(E, ..., Yname=Yname))
}
  
