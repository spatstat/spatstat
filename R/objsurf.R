#
#  objsurf.R
#
#  surface of the objective function for an M-estimator
#
#  $Revision: 1.4 $ $Date: 2015/10/21 09:06:57 $
#

objsurf <- function(x, ...) {
  UseMethod("objsurf")
}

objsurf.kppm <- objsurf.dppm <- function(x, ..., ngrid=32, ratio=1.5, verbose=TRUE) {
  Fit <- x$Fit
  switch(Fit$method,
         mincon = {
           result <- objsurf(Fit$mcfit, ...,
                             ngrid=ngrid, ratio=ratio, verbose=verbose)
         },
         clik = {
           optpar  <- x$par
           objfun  <- Fit$objfun
           objargs <- Fit$objargs
           result  <- objsurfEngine(objfun, optpar, objargs, ...,
                                    ngrid=ngrid, ratio=ratio, verbose=verbose)
         })
  return(result)
}

objsurf.minconfit <- function(x, ..., ngrid=32, ratio=1.5, verbose=TRUE) {
  optpar  <- x$par.canon %orifnull% x$par
  objfun  <- x$objfun
  objargs <- x$objargs
  dotargs <- x$dotargs
  objsurfEngine(objfun, optpar, objargs, ...,
                dotargs=dotargs,
                ngrid=ngrid, ratio=ratio, verbose=verbose)
}

objsurfEngine <- function(objfun, optpar, objargs, 
                          ...,
                          dotargs=list(),
                          objname="objective", 
                          ngrid=32, ratio=1.5, verbose=TRUE) {
  trap.extra.arguments(...)
  if(!is.function(objfun))
    stop("Object is in an outdated format and needs to be re-fitted")
  npar    <- length(optpar)
  if(npar != 2)
    stop("Only implemented for functions of 2 arguments")
  # create grid of parameter values
  ratio <- ensure2vector(ratio)
  ngrid <- ensure2vector(ngrid)
  stopifnot(all(ratio > 1))
  xgrid <- seq(optpar[1]/ratio[1], optpar[1] * ratio[1], length=ngrid[1])
  ygrid <- seq(optpar[2]/ratio[2], optpar[2] * ratio[2], length=ngrid[2])
  pargrid <- expand.grid(xgrid, ygrid)
  colnames(pargrid) <- names(optpar)
  # evaluate
  if(verbose) cat(paste("Evaluating", nrow(pargrid), "function values..."))
  values <- do.call("apply",
                    append(list(pargrid, 1, objfun, objargs=objargs), dotargs))
  if(verbose) cat("Done.\n")
  result <- list(x=xgrid, y=ygrid, z=matrix(values, ngrid[1], ngrid[2]))
  attr(result, "optpar") <- optpar
  attr(result, "objname") <- "contrast"
  class(result) <- "objsurf"
  return(result)
}

print.objsurf <- function(x, ...) {
  cat("Objective function surface\n")
  optpar <- attr(x, "optpar")
  objname <- attr(x, "objname")
  nama <- names(optpar)
  cat("Parameter ranges:\n")
  cat(paste(paste0(nama[1], ":"), prange(range(x$x)), "\n"))
  cat(paste(paste0(nama[2], ":"), prange(range(x$y)), "\n"))
  cat(paste("Function value:", objname, "\n"))
  invisible(NULL)
}

image.objsurf <- plot.objsurf <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  optpar <- attr(x, "optpar")
  nama <- names(optpar)
  do.call("image",
          resolve.defaults(list(x=unclass(x)),
                           list(...),
                           list(xlab=nama[1], ylab=nama[2], main=xname)))
  abline(v=optpar[1], lty=3)
  abline(h=optpar[2], lty=3)
  invisible(NULL)
}

contour.objsurf <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  optpar <- attr(x, "optpar")
  nama <- names(optpar)
  do.call("contour",
          resolve.defaults(list(x=unclass(x)),
                           list(...),
                           list(xlab=nama[1], ylab=nama[2], main=xname)))
  abline(v=optpar[1], lty=3)
  abline(h=optpar[2], lty=3)
  invisible(NULL)
}

  
persp.objsurf <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  optpar <- attr(x, "optpar")
  objname <- attr(x, "objname")
  nama <- names(optpar)
  r <- do.call("persp",
               resolve.defaults(list(x=x$x, y=x$y, z=x$z),
                                list(...),
                                list(xlab=nama[1], ylab=nama[2],
                                     zlab=objname, main=xname)))
  invisible(r)
}


