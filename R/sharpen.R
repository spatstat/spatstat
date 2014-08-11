#
#      sharpen.R
#
#      $Revision: 1.5 $  $Date: 2010/11/27 01:52:33 $
#

sharpen <- function(X, ...) {
  UseMethod("sharpen")
}

sharpen.ppp <- function(X, sigma=NULL, ..., varcov=NULL,
                        edgecorrect=FALSE) {
  stopifnot(is.ppp(X))
  Yx <- smooth.ppp(X %mark% X$x,
                   at="points", sigma=sigma, varcov=varcov, edge=TRUE)
  Yy <- smooth.ppp(X %mark% X$y,
                   at="points", sigma=sigma, varcov=varcov, edge=TRUE)
  # trap NaN etc
  nbad <- sum(!(is.finite(Yx) & is.finite(Yy)))
  if(nbad > 0)
    stop(paste(nbad,
               ngettext(nbad, "point is", "points are"),
               "undefined due to numerical problems;",
               "smoothing parameter is probably too small"))
  #
  W <- as.owin(X)
  if(edgecorrect) {
    # convolve x and y coordinate functions with kernel
    xim <- as.im(function(x,y){x}, W)
    yim <- as.im(function(x,y){y}, W)
    xblur <- blur(xim, sigma=sigma, varcov=varcov, normalise=TRUE, ...)
    yblur <- blur(yim, sigma=sigma, varcov=varcov, normalise=TRUE, ...)
    # evaluate at data locations 
    xx <- safelookup(xblur, X, warn=FALSE)
    yy <- safelookup(yblur, X, warn=FALSE)
    # estimated vector bias of sharpening procedure
    xbias <- xx - X$x
    ybias <- yy - X$y
    # adjust
    Yx <- Yx - xbias
    Yy <- Yy - ybias
    # check this does not place points outside window
    if(any(uhoh <- !inside.owin(Yx, Yy, W))) {
      # determine mass of edge effect
      edgeim <- blur(as.im(W), sigma=sigma, varcov=varcov, normalise=FALSE, ...)
      edg <- safelookup(edgeim, X[uhoh], warn=FALSE)
      # contract bias correction
      Yx[uhoh] <- (1 - edg) * X$x[uhoh] + edg * Yx[uhoh]
      Yy[uhoh] <- (1 - edg) * X$y[uhoh] + edg * Yy[uhoh]
    }
    # check again
    if(any(nbg <- !inside.owin(Yx, Yy, W))) {
      # give up
      Yx[nbg] <- X$x[nbg]
      Yy[nbg] <- X$y[nbg]
    }
  }
  # make point pattern
  Y <- ppp(Yx, Yy, marks=marks(X), window=W)
  # tack on smoothing information
  attr(Y, "sigma") <- sigma
  attr(Y, "varcov") <- varcov
  attr(Y, "edgecorrected") <- edgecorrect
  return(Y)
}
