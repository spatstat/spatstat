#
# bw.optim.R
#
#  Class of optimised bandwidths
#  Plotting the object displays the optimisation criterion
#
#  $Revision: 1.4 $  $Date: 2011/09/12 08:22:08 $
#

bw.optim <- function(cv, h, iopt=which.min(cv), ...,
                     xlab="h", ylab="CV") {
  stopifnot(is.numeric(cv))
  stopifnot(is.numeric(h))
  stopifnot(length(h) == length(cv))
  result <- h[iopt]
  attr(result, "cv") <- cv
  attr(result, "h") <- h
  attr(result, "iopt") <- iopt
  attr(result, "labels") <- list(xlab=xlab,
                                 ylab=ylab)
  attr(result, "info") <- list(...)
  class(result) <- "bw.optim"
  return(result)
}

print.bw.optim <- function(x, ...) {
  y <- as.numeric(x)
  names(y) <- attr(x, "labels")$xlab
  print(y, ...)
  return(invisible(NULL))
}

plot.bw.optim <- function(x, ..., add=FALSE) {
  xname <- short.deparse(substitute(x))
  h <- attr(x, "h")
  cv <- attr(x, "cv")
  iopt <- attr(x, "iopt")
  labels <- attr(x, "labels")
  if(add) lines(x=h, y=cv, ...) else  
    do.call("plot",
            resolve.defaults(list(x=h, y=cv),
                             list(...),
                             list(main=xname, type="l"),
                             labels))
  return(invisible(NULL))
}


