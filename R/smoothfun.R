##
## smoothfun.R
##
## Exact 'funxy' counterpart of Smooth.ppp
##
##  $Revision: 1.6 $ $Date: 2018/04/13 06:24:43 $


Smoothfun <- function(X, ...) {
  UseMethod("Smoothfun")
}

Smoothfun.ppp <- function(X, sigma=NULL, ...,
                          weights=NULL, edge=TRUE, diggle=FALSE) {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE))
    stop("X should be a marked point pattern")
  ## handle weights now
  weightsgiven <- !missing(weights) && !is.null(weights) 
  if(weightsgiven) {
    # convert to numeric
    if(is.im(weights)) {
      weights <- safelookup(weights, X) # includes warning if NA
    } else if(is.expression(weights)) 
      weights <- eval(weights, envir=as.data.frame(X), enclos=parent.frame())
    if(length(weights) == 0)
      weightsgiven <- FALSE
  }
  if(weightsgiven) {
    check.nvector(weights, npoints(X))
  } else weights <- NULL
  ## 
  X <- coerce.marks.numeric(X)
  ## 
  stuff <- list(X=X, weights=weights, edge=edge, diggle=diggle)
  ## 
  ## determine smoothing parameters
  ker <- resolve.2D.kernel(sigma=sigma, ...,
                           x=X, bwfun=bw.smoothppp, allow.zero=TRUE)
  stuff <- append(stuff, ker[c("sigma", "varcov")])
  ##
  g <- function(x, y=NULL) {
    Y <- xy.coords(x, y)[c("x", "y")]
    with(stuff,
         smoothcrossEngine(Xdata=X,
                           Xquery=as.ppp(Y, X$window),
                           values=marks(X),
                           sigma=sigma,
                           varcov=varcov, 
                           weights=weights,
                           edge=edge, diggle=diggle))
  }
  g <- funxy(g, as.rectangle(as.owin(X)))
  class(g) <- c("Smoothfun", class(g))
  return(g)
}

print.Smoothfun <- function(x, ...) {
  cat("function(x,y)", "which returns",
      "values", "interpolated from", fill=TRUE)
  X <- get("X", envir=environment(x))
  print(X, ...)
  return(invisible(NULL))
}

## Method for as.im
## (enables plot.funxy, persp.funxy, contour.funxy to work for this class)

as.im.Smoothfun <- function(X, W=Window(X), ..., approx=TRUE) {
  stuff <- get("stuff", envir=environment(X))
  if(!approx) {
    #' evaluate exactly at grid points 
    result <- as.im.function(X, W=W, ...)
  } else {
    #' faster, approximate evaluation using FFT
    if(!is.null(W)) stuff$X <- stuff$X[W]
    result <- do.call(Smooth, resolve.defaults(list(...), stuff))
  }
  return(result)
}


  
