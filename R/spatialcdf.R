##
## spatialcdf.R
##
##  $Revision: 1.1 $ $Date: 2014/06/26 03:11:49 $
##

spatialcdf <- function(Z, weights=NULL, normalise=FALSE, ...,
                       W=NULL, Zname=NULL) {
  Zdefaultname <- singlestring(short.deparse(substitute(Z)))
  if(is.character(Z) && length(Z) == 1) {
    if(is.null(Zname)) Zname <- Z
    switch(Zname,
           x={
             Z <- function(x,y) { x }
           }, 
           y={
             Z <- function(x,y) { y }
           },
           stop("Unrecognised covariate name")
         )
  }
  if(is.null(Zname)) Zname <- Zdefaultname
  ##
  if(is.ppm(weights) || is.kppm(weights) || is.dppm(weights)) {
    Q <- quad.ppm(as.ppm(weights))
    loc <- as.ppp(Q)
    df <- mpl.get.covariates(list(Z=Z), loc, covfunargs=list(...))
    df$wt <- fitted(weights) * w.quad(Q)
    wtname <- if(normalise) "fraction of points" else "number of points"
  } else {
    if(is.null(W)) W <- as.owin(weights, fatal=FALSE)
    if(is.null(W)) W <- as.owin(Z, fatal=FALSE)
    if(is.null(W)) stop("No information specifying the spatial window")
    if(is.null(weights)) weights <- 1
    M <- as.mask(W, ...)
    loc <- rasterxy.mask(M, drop=TRUE)
    df <- mpl.get.covariates(list(Z=Z, weights=weights), loc,
                             covfunargs=list(...))
    pixelarea <- with(unclass(M), xstep * ystep)
    df$wt <- rep(pixelarea, nrow(df))
    wtname <- if(normalise) "fraction of weight" else "weight"
  }
  if(normalise) 
    df$wt <- with(df, wt/sum(wt))
  G <- with(df, ewcdf(Z, wt))
  class(G) <- c("spatialcdf", class(G))
  attr(G, "call") <- sys.call()
  attr(G, "Zname") <- Zname
  attr(G, "ylab") <- paste("Cumulative", wtname)
  return(G)
}

plot.spatialcdf <- function(x, ..., xlab, ylab) {
  if(missing(xlab) || is.null(xlab))
    xlab <- attr(x, "Zname")
  if(missing(ylab) || is.null(ylab))
    ylab <- attr(x, "ylab")
  plot.ecdf(x, ..., xlab=xlab, ylab=ylab)
}

