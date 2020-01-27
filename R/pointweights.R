#'  pointweights.R
#'
#'  get a valid vector of weights for a point pattern
#'
#'  Argument 'weights' is usually passed from a user-level function
#'  It may be:
#'        a numeric vector
#'        a single number
#'        a function(x,y)
#'        a pixel image
#'        an expression involving the coordinates and marks
#' 
#'  $Revision: 1.2 $ $Date: 2020/01/27 09:08:06 $

pointweights <- function(X, ..., weights=NULL, parent=NULL) {
  if(is.null(weights)) return(NULL)
  nX <- npoints(X)
  if(is.numeric(weights) && is.vector(as.numeric(weights))) {
    if(length(weights) == 1) weights <- rep(weights, nX)
  } else if(is.im(weights)) {
    weights <- safelookup(weights, X) # includes warning if NA
  } else if(is.function(weights)) {
      weights <- weights(X$x, X$y)
  } else if(is.expression(weights)) {
    #' evaluate expression in data frame of coordinates and marks
    df <- as.data.frame(X)
    weights <- try(eval(weights, envir=df, enclos=parent))
    if(inherits(weights, "try-error"))
      stop("Unable to evaluate expression for weights", call.=FALSE)
    if(length(weights) == 0) return(NULL)
  } else stop(paste("Argument 'weights' should be",
                    "a numeric vector, a function, an image,",
                    "or an expression"), call.=FALSE)
  check.nvector(weights, nX)
  return(weights)
}

