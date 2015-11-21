#'
#'   zclustermodel.R
#'
#' Experimental
#'

zclustermodel <- function(name="Thomas", ..., mu, kappa, scale) {
  if(missing(kappa)) stop("The parent intensity kappa must be given")
  if(missing(mu)) stop("The mean cluster size mu must be given")
  if(missing(scale)) stop("The cluster scale must be given")
  rules <- spatstatClusterModelInfo(name)
  argh <- list(startpar=c(kappa=kappa, scale=scale), ...)
  argh <- do.call(rules$resolvedots, argh)
  par <- argh$startpar
  other <- argh[names(argh) != "startpar"]
  clustargs <- rules$checkclustargs(other$margs, old=FALSE)
  out <- list(name=name, rules=rules, par=par, mu=mu, clustargs=clustargs,
              other=other)
  class(out) <- "zclustermodel"
  return(out)
}

print.zclustermodel <- local({

  print.zclustermodel <- function(x, ...) {
    with(x, {
      splat(rules$printmodelname(list(clustargs=clustargs)))
      newpar <- rules$checkpar(par, old=FALSE)
      splat("Parent intensity kappa =", blurb("kappa", newpar["kappa"]))
      splat("Cluster scale = ", newpar["scale"])
      splat("Mean cluster size mu =", blurb("mu", mu))
      if(length(clustargs) > 0) {
        hdr <- paste("Cluster shape",
                     ngettext(length(clustargs), "parameter:", "parameters:"))
        if(is.list(clustargs) &&
           all(sapply(clustargs, is.numeric)) &&
           all(sapply(clustargs, length) == 1)) {
          splat(hdr,
                paste(names(clustargs), as.numeric(clustargs),
                      sep="=",
                      collapse=", "))
        } else {
          splat(hdr)
          print(clustargs)
        }
      }
    })
    return(invisible(NULL))
  }

  blurb <- function(name, value) {
    if(is.numeric(value)) as.character(value) else
    if(is.im(value)) "[image]" else "[unrecognized format]"
  }

  print.zclustermodel
})

                             
pcfmodel.zclustermodel <- function(model, ...) {
  p <- model$rules$pcf
  mpar <- model$par
  other <- model$other
  f <- function(r) {
    do.call(p, c(list(par=mpar, rvals=r), other, model$rules["funaux"]))
  }
  return(f)
}

predict.zclustermodel <- function(object, ...,
                                 locations,
                                 type="intensity",
                                 ngrid=NULL) {
  ## limited use!!!
  if(!identical(type, "intensity"))
    stop("Sorry, only type='intensity' is implemented")
  lambda <- object$par["kappa"] * object$mu
  if(is.numeric(lambda)) {
    if(is.ppp(locations))
      return(rep(lambda, npoints(locations)))
    W <- as.owin(locations)
    if(!is.mask(W))
      W <- as.mask(W, dimyx=ngrid, ...)
    return(as.im(lambda, W=W))
  }
  return(lambda[locations])
}


