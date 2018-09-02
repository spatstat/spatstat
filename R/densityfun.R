##
## densityfun.R
##
## Exact 'funxy' counterpart of density.ppp
##
##  $Revision: 1.8 $ $Date: 2018/09/02 07:47:31 $


densityfun <- function(X, ...) {
  UseMethod("densityfun")
}

densityfun.ppp <- function(X, sigma=NULL, ...,
                          weights=NULL, edge=TRUE, diggle=FALSE) {
  verifyclass(X, "ppp")
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
  stuff <- list(Xdata=X, weights=weights, edge=edge, diggle=diggle, ...)
  ##
  ## determine smoothing parameters
  ker <- resolve.2D.kernel(sigma=sigma, ...,
                           x=X, bwfun=bw.diggle, allow.zero=TRUE)
  stuff[c("sigma", "varcov")]  <- ker[c("sigma", "varcov")]
  ##
  g <- function(x, y=NULL) {
    Y <- xy.coords(x, y)[c("x", "y")]
    Xquery <- as.ppp(Y, Window(stuff$Xdata))
    do.call(densitycrossEngine, append(list(Xquery=Xquery), stuff))
  }
  g <- funxy(g, as.rectangle(as.owin(X)))
  class(g) <- c("densityfun", class(g))
  return(g)
}

print.densityfun <- function(x, ...) {
  cat("function(x,y)", "which returns",
      "kernel estimate of intensity for", fill=TRUE)
  X <- get("X", envir=environment(x))
  print(X, ...)
  return(invisible(NULL))
}

## Method for as.im
## (enables plot.funxy, persp.funxy, contour.funxy to work for this class)

as.im.densityfun <- function(X, W=Window(X), ..., approx=TRUE) {
  if(!approx) {
    #' evaluate exactly at grid points using as.im.funxy -> as.im.function
    result <- as.im.function(X, W=W, ...)
  } else {
    #' faster, approximate evaluation using FFT
    stuff <- get("stuff", envir=environment(X))
    if(!missing(W)) stuff$X <- stuff$X[W]
    names(stuff)[names(stuff) == "X"] <- "x"
    result <- do.call(density, resolve.defaults(list(...), stuff))
  }
  return(result)
}


  
