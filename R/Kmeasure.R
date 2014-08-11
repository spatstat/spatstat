#
#           Kmeasure.R
#
#           $Revision: 1.42 $    $Date: 2013/04/25 06:37:43 $
#
#     Kmeasure()         compute an estimate of the second order moment measure
#
#     Kest.fft()        use Kmeasure() to form an estimate of the K-function
#
#     second.moment.calc()    underlying algorithm
#

Kmeasure <- function(X, sigma, edge=TRUE, ..., varcov=NULL) {
  stopifnot(is.ppp(X))
  
  sigma.given <- !missing(sigma) && !is.null(sigma)
  varcov.given <- !is.null(varcov)
  ngiven <- sigma.given + varcov.given
  if(ngiven == 2)
    stop(paste("Give only one of the arguments",
               sQuote("sigma"), "and", sQuote("varcov")))
  if(ngiven == 0)
    stop(paste("Please specify smoothing bandwidth", sQuote("sigma"),
               "or", sQuote("varcov")))
  if(varcov.given) {
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov)==2 )
    sigma <- NULL
  } else {
    stopifnot(is.numeric(sigma))
    stopifnot(length(sigma) %in% c(1,2))
    stopifnot(all(sigma > 0))
    if(length(sigma) == 2) {
      varcov <- diag(sigma^2)
      sigma <- NULL
    }
  }  

  second.moment.calc(X, sigma=sigma, edge, "Kmeasure", varcov=varcov)
}

second.moment.calc <- function(x, sigma=NULL, edge=TRUE,
                               what="Kmeasure", debug=FALSE, ...,
                               varcov=NULL, expand=FALSE) {
  if(is.null(sigma) && is.null(varcov))
    stop("must specify sigma or varcov")
  choices <- c("kernel", "smooth", "Kmeasure", "Bartlett", "edge", "all", "smoothedge")
  if(!(what %in% choices))
    stop(paste("Unknown choice: what = ", sQuote(what),
               "; available options are:",
               paste(sQuote(choices), collapse=", ")))
  sig <- if(!is.null(sigma)) sigma else max(c(diag(varcov), sqrt(det(varcov))))
  win <- rescue.rectangle(as.owin(x))
  rec <- as.rectangle(win)
  across <- min(diff(rec$xrange), diff(rec$yrange))
  # determine whether to expand window
  if(!expand || (6 * sig < across))
    result <- second.moment.engine(x, sigma=sigma, edge=edge,
                                   what=what, debug=debug, ..., varcov=varcov)
  else {
    # need to expand window
    bigger <- grow.rectangle(rec, (7 * sig - across)/2)
    if(is.ppp(x)) {
      # pixellate first (to preserve pixel resolution)
      X <- pixellate(x, ..., padzero=TRUE)
      np <- npoints(x)
    } else if(is.im(x)) {
      X <- x
      np <- NULL
    } else stop("x should be an image or a point pattern")
    # now expand 
    X <- rebound.im(X, bigger)
    X <- na.handle.im(X, 0)
    out <- second.moment.engine(X, sigma=sigma, edge=edge,
                                what=what, debug=debug, ...,
                                obswin=win, varcov=varcov, npts=np)
    # now clip it
    fbox <- shift(rec, origin="midpoint")
    result <- switch(what,
                     kernel   = out[fbox],
                     smooth   = out[win],
                     Kmeasure = out[fbox],
                     Bartlett = out[fbox],
                     edge     = out[win],
                     smoothedge = list(smooth=out$smooth[win],
                                       edge  =out$edge[win]),
                     all      =
                     list(kernel=out$kernel[fbox],
                          smooth=out$smooth[win],
                          Kmeasure=out$Kmeasure[fbox],
                          Bartlett=out$Bartlett[fbox],
                          edge=out$edge[win]))
  }
  return(result)
}

second.moment.engine <- function(x, sigma=NULL, edge=TRUE,
                                 what="Kmeasure", debug=FALSE, ...,
                                 obswin = as.owin(x), varcov=NULL,
                                 npts=NULL)
{
  stopifnot(what %in% c("kernel", "smooth", "Kmeasure",
                        "Bartlett", "edge", "all", "smoothedge"))
  if(is.ppp(x)) 
    # convert list of points to mass distribution
    X <- pixellate(x, ..., padzero=TRUE)
  else if(is.im(x))
    X <- x
  else stop("internal error: unrecognised format for x")
  # ensure obswin has same bounding frame as X
  if(!missing(obswin))
    obswin <- rebound.owin(obswin, as.rectangle(X))
  #
  if(is.null(npts) && is.ppp(x))
      npts <- npoints(x)
  # go to work
  Y <- X$v
  xw <- X$xrange
  yw <- X$yrange
  # pad with zeroes
  nr <- nrow(Y)
  nc <- ncol(Y)
  Ypad <- matrix(0, ncol=2*nc, nrow=2*nr)
  Ypad[1:nr, 1:nc] <- Y
  lengthYpad <- 4 * nc * nr
  # corresponding coordinates
  xw.pad <- xw[1] + 2 * c(0, diff(xw))
  yw.pad <- yw[1] + 2 * c(0, diff(yw))
  xcol.pad <- X$xcol[1] + X$xstep * (0:(2*nc-1))
  yrow.pad <- X$yrow[1] + X$ystep * (0:(2*nr-1))
  # set up Gauss kernel
  xcol.G <- X$xstep * c(0:(nc-1),-(nc:1))
  yrow.G <- X$ystep * c(0:(nr-1),-(nr:1))
  xx <- matrix(xcol.G[col(Ypad)], ncol=2*nc, nrow=2*nr)
  yy <- matrix(yrow.G[row(Ypad)], ncol=2*nc, nrow=2*nr)
  if(!is.null(sigma)) {
#    if(max(abs(diff(xw)),abs(diff(yw))) < 6 * sigma)
#      warning("sigma is too large for this window")
    Kern <- exp(-(xx^2 + yy^2)/(2 * sigma^2))/(2 * pi * sigma^2) * X$xstep * X$ystep
  } else if(!is.null(varcov)) {
    # anisotropic kernel
    detSigma <- det(varcov)
    Sinv <- solve(varcov)
    const <- X$xstep * X$ystep/(2 * pi * sqrt(detSigma))
    Kern <- const * exp(-(xx * (xx * Sinv[1,1] + yy * Sinv[1,2])
                          + yy * (xx * Sinv[2,1] + yy * Sinv[2,2]))/2)
  } else 
    stop("Must specify either sigma or varcov")

  # these options call for several image outputs
  if(what %in% c("all", "smoothedge"))
    result <- list()
  
  if(what %in% c("kernel", "all")) {
    # return the kernel
    # first rearrange it into spatially sensible order (monotone x and y)
    rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
    ctwist <- (-nc):(nc-1) %% (2*nc) + 1
    if(debug) {
      if(any(fave.order(xcol.G) != rtwist))
        cat("something round the twist\n")
    }
    Kermit <- Kern[ rtwist, ctwist]
    ker <- im(Kermit, xcol.G[ctwist], yrow.G[ rtwist], unitname=unitname(x))
    if(what == "kernel")
      return(ker)
    else 
      result$kernel <- ker
  }
  # convolve using fft
  fK <- fft(Kern)
  if(what != "edge") {
    fY <- fft(Ypad)
    sm <- fft(fY * fK, inverse=TRUE)/lengthYpad
    if(debug) {
      cat(paste("smooth: maximum imaginary part=",
                signif(max(Im(sm)),3), "\n"))
      if(!is.null(npts))
        cat(paste("smooth: mass error=",
                  signif(sum(Mod(sm))-npts,3), "\n"))
    }
  }
  if(what %in% c("smooth", "all", "smoothedge")) {
    # return the smoothed point pattern without edge correction
    smo <- im(Re(sm)[1:nr, 1:nc], xcol.pad[1:nc], yrow.pad[1:nr],
              unitname=unitname(x))
    if(what == "smooth")
      return(smo)
    else
      result$smooth <- smo
  }

  if(what != "edge")
    bart <- Mod(fY)^2 * fK
  
  if(what %in% c("Bartlett", "all")) {
     # return Bartlett spectrum
     # rearrange into spatially sensible order (monotone x and y)
    rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
    ctwist <- (-nc):(nc-1) %% (2*nc) + 1
    Bart <- bart[ rtwist, ctwist]
    Bartlett <- im(Mod(Bart),(-nc):(nc-1), (-nr):(nr-1))
    if(what == "Bartlett")
      return(Bartlett)
    else
      result$Bartlett <- Bartlett
  }
  
  #### ------- Second moment measure --------------
  #
  if(what != "edge") {
    mom <- fft(bart, inverse=TRUE)/lengthYpad
    if(debug) {
      cat(paste("2nd moment measure: maximum imaginary part=",
                signif(max(Im(mom)),3), "\n"))
      if(!is.null(npts))
        cat(paste("2nd moment measure: mass error=",
                  signif(sum(Mod(mom))-npts^2, 3), "\n"))
    }
    mom <- Mod(mom)
    # subtract (delta_0 * kernel) * npts
    if(is.null(npts))
      stop("Internal error: second moment measure requires npts")
    mom <- mom - npts* Kern
  }
  # edge correction
  if(edge || what %in% c("edge", "all", "smoothedge")) {
    # compute kernel-smoothed set covariance
    M <- as.mask(obswin, dimyx=c(nr, nc))$m
    # previous line ensures M has same dimensions and scale as Y 
    Mpad <- matrix(0, ncol=2*nc, nrow=2*nr)
    Mpad[1:nr, 1:nc] <- M
    lengthMpad <- 4 * nc * nr
    fM <- fft(Mpad)
    if(edge && what != "edge") {
      co <- fft(Mod(fM)^2 * fK, inverse=TRUE)/lengthMpad
      co <- Mod(co) 
      a <- sum(M)
      wt <- a/co
      me <- spatstat.options("maxedgewt")
      weight <- matrix(pmin.int(me, wt), ncol=2*nc, nrow=2*nr)
#      if(debug) browser()
      #
      # apply edge correction
      mom <- mom * weight
      # set to NA outside 'reasonable' region
      mom[wt > 10] <- NA
    }
  }
  if(what != "edge") {
    # rearrange into spatially sensible order (monotone x and y)
    rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
    ctwist <- (-nc):(nc-1) %% (2*nc) + 1
    mom <- mom[ rtwist, ctwist]
    if(debug) {
      if(any(fave.order(xcol.G) != rtwist))
        cat("something round the twist\n")
    }
  }
  if(what %in% c("edge", "all", "smoothedge")) {
    # return convolution of window with kernel
    # (evaluated inside window only)
    con <- fft(fM * fK, inverse=TRUE)/lengthMpad
    edg <- Mod(con[1:nr, 1:nc])
    edg <- im(edg, xcol.pad[1:nc], yrow.pad[1:nr], unitname=unitname(x))
    if(what == "edge") 
      return(edg)
    else
      result$edge <- edg
  }
  if(what == "smoothedge")
    return(result)
  # Second moment measure, density estimate
  # Divide by number of points * lambda and convert mass to density
  pixarea <- with(X, xstep * ystep)
  mom <- mom * area.owin(obswin) / (pixarea * npts^2)
  # this is the second moment measure
  mm <- im(mom, xcol.G[ctwist], yrow.G[rtwist], unitname=unitname(x))
  if(what == "Kmeasure")
    return(mm)
  else 
    result$Kmeasure <- mm
  # what = "all", so return all computed objects
  return(result)
}

Kest.fft <- function(X, sigma, r=NULL, breaks=NULL) {
  verifyclass(X, "ppp")
  W <- X$window
  lambda <- X$n/area.owin(W)
  rmaxdefault <- rmax.rule("K", W, lambda)        
  bk <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  breaks <- bk$val
  rvalues <- bk$r
  u <- Kmeasure(X, sigma)
  xx <- rasterx.im(u)
  yy <- rastery.im(u)
  rr <- sqrt(xx^2 + yy^2)
  tr <- whist(rr, breaks, u$v)
  K  <- cumsum(tr)
  rmax <- min(rr[is.na(u$v)])
  K[rvalues >= rmax] <- NA
  result <- data.frame(r=rvalues, theo=pi * rvalues^2, border=K)
  w <- X$window
  alim <- c(0, min(diff(w$xrange), diff(w$yrange))/4)
  out <- fv(result,
            "r", substitute(K[fft](r), NULL),
            "border",
             . ~ r, alim,
            c("r", "{%s^{pois}}(r)", "{hat(%s)^{bord}}(r)"),
            c("distance argument r",
              "theoretical Poisson %s",
              "border-corrected estimate of %s"),
            fname="K[fft]",
            unitname=unitname(X)
            )
  return(out)
}

