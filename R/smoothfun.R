##
## smoothfun.R
##
## Exact 'funxy' counterpart of Smooth.ppp
##
##  $Revision: 1.2 $ $Date: 2016/02/11 10:17:12 $


Smoothfun <- function(X, ...) {
  UseMethod("Smoothfun")
}

Smoothfun.ppp <- function(X, sigma=NULL, ...,
                          weights=NULL, edge=TRUE, diggle=FALSE) {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE))
    stop("X should be a marked point pattern")
  stuff <- list(X=X, weights=weights, edge=edge, diggle=diggle)
  X <- coerce.marks.numeric(X)
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

as.im.Smoothfun <- function(X, W=NULL, ...) {
  stuff <- get("stuff", envir=environment(X))
  if(!is.null(W)) stuff$X <- stuff$X[W]
  do.call(Smooth, resolve.defaults(list(...), stuff))
}


  
