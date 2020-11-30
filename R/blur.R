#
# blur.R
#
# apply Gaussian blur to an image
#
#    $Revision: 1.25 $   $Date: 2020/11/30 07:16:06 $
#

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
  
