#
# blur.R
#
# apply Gaussian blur to an image
#
#    $Revision: 1.23 $   $Date: 2020/02/04 07:05:28 $
#
fillNA <- function(x, value=0) {
  stopifnot(is.im(x))
  v <- x$v
  v[is.na(v)] <- value
  x$v <- v
  return(x)
}

Smooth.im <- function(X, sigma=NULL, ..., kernel="gaussian",
                      normalise=FALSE, bleed=TRUE, varcov=NULL) {
  blur(X, sigma=sigma, ..., kernel=kernel,
       normalise=normalise, bleed=bleed, varcov=varcov)
}

blur <- function(x, sigma=NULL, ..., kernel="gaussian",
                 normalise=FALSE, bleed=TRUE, varcov=NULL) {
  stopifnot(is.im(x))
  # determine smoothing kernel 
  sigma.given <- !is.null(sigma)
  varcov.given <- !is.null(varcov)
  if (sigma.given) {
    stopifnot(is.numeric(sigma))
    stopifnot(length(sigma) %in% c(1, 2))
    stopifnot(all(sigma > 0))
  }
  if (varcov.given)
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov) == 2)
  ngiven <- varcov.given + sigma.given
  switch(ngiven + 1L,
         {
           sigma <- (1/8) * min(diff(x$xrange), diff(x$yrange))
         }, {
           if (sigma.given && length(sigma) == 2)
             varcov <- diag(sigma^2)
           if (!is.null(varcov))
             sigma <- NULL
         }, {
           stop(paste("Give only one of the arguments", sQuote("sigma"),
                      "and", sQuote("varcov")))
         })
  # replace NA's in image raster by zeroes 
  X <- fillNA(x, 0)
  # convolve with Gaussian
  Y <- second.moment.calc(X, sigma=sigma, ..., kernel=kernel,
                          varcov=varcov, what="smooth")
  # if no bleeding, we restrict data to the original boundary
  if(!bleed)
    Y$v[is.na(x$v)] <- NA
  # 
  if(!normalise)
    return(Y)
  # normalisation:
  # convert original image to window (0/1 image)
  Xone <- x
  isna <- is.na(x$v)
  Xone$v[isna] <- 0
  Xone$v[!isna] <- 1
  # convolve with Gaussian
  Ydenom <- second.moment.calc(Xone, sigma=sigma, ..., kernel=kernel,
                               varcov=varcov, what="smooth")
  # normalise
  Z <- eval.im(Y/Ydenom)
  return(Z)
}
  
safelookup <- function(Z, x, factor=2, warn=TRUE) {
  #' x is a ppp
  #' evaluates Z[x], replacing any NA's by blur(Z)[x]
  Zvals <- Z[x, drop=FALSE]
  isna <- is.na(Zvals)
  if(!any(isna))
     return(Zvals)
  #' First pass - look up values at neighbouring pixels if valid
  Xbad <- x[isna]
  rc <- nearest.valid.pixel(Xbad$x, Xbad$y, Z)
  Nvals <- Z$v[cbind(rc$row, rc$col)]
  fixed <- !is.na(Nvals)
  Zvals[isna] <- Nvals
  if(all(fixed))
    return(Zvals)
  isna[isna] <- !fixed
  #' Second pass 
  Xbad <- x[isna]
  #' expand domain of Z 
  RX <- as.rectangle(x)
  RZ <- as.rectangle(Z)
  bb <- boundingbox(RX, RZ)
  pixdiam <- sqrt(Z$xstep^2 + Z$ystep^2)
  big <- grow.rectangle(bb, 2 * pixdiam)
  Z <- rebound.im(Z, big)
  #'
  isfac <- (Z$type == "factor")
  if(!isfac) {
    #' Numerical extrapolation: blur by a few pixels
    Zblur <- blur(Z, factor * pixdiam, bleed=TRUE, normalise=TRUE)
    Bvals <- Zblur[Xbad, drop=FALSE]
    Zvals[isna] <- Bvals
    fixed <- !is.na(Bvals)
    if(warn && any(fixed))
      warning(paste("Values for", sum(fixed), 
                    "points lying slightly outside the pixel image domain",
                    "were estimated by convolution"),
              call.=FALSE)
    if(all(fixed))
      return(Zvals)
    isna[isna] <- notfixed <- !fixed
    Xbad <- Xbad[notfixed]
  }
  #' Third pass
  #' last resort: project to nearest pixel at any distance
  W <- as.mask(Z)
  eW <- exactPdt(W)
  ## discretise points of Xbad
  Gbad <- nearest.raster.point(Xbad$x, Xbad$y, W)
  ijGbad <- cbind(Gbad$row, Gbad$col)
  ## find nearest pixels inside domain
  iclosest <- eW$row[ijGbad]
  jclosest <- eW$col[ijGbad]
  ## look up values of Z
  Cvals <- Z$v[cbind(iclosest, jclosest)]
  fixed <- !is.na(Cvals)
  Zvals[isna] <- Cvals
  if(warn && any(fixed)) 
    warning(paste(if(isfac) "Categorical values" else "Values",
                  "for", sum(fixed), "points lying",
                  if(isfac) "outside" else "far outside",
                  "the pixel image domain",
                  "were estimated by projection to the nearest pixel"),
            call.=FALSE)
  if(!all(fixed))
    stop(paste("Internal error:", sum(!fixed),
               "pixel values were NA, even after projection"),
         call.=FALSE)
  return(Zvals)
}

nearestValue <- function(X) {
  #' for each raster location, look up the nearest defined pixel value
  X <- as.im(X)
  if(!anyNA(X)) return(X)
  Y <- X ## copy dimensions, value type, units etc etc
  W <- as.mask(X)
  eW <- exactPdt(W)
  iclosest <- as.vector(eW$row)
  jclosest <- as.vector(eW$col)
  ## look up values of Z
  Y$v[] <- X$v[cbind(iclosest, jclosest)]
  return(Y)
}

