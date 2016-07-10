#
#           Kmeasure.R
#
#           $Revision: 1.54 $    $Date: 2016/07/10 08:06:15 $
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

  second.moment.calc(x=X, sigma=sigma, edge=edge,
                     what="Kmeasure", varcov=varcov, ...)
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

  xtype <- if(is.ppp(x)) "ppp" else
           if(is.im(x)) "im" else 
           if(all(unlist(lapply(x, is.im)))) "imlist" else
           stop("x should be a point pattern or a pixel image")

  nimages <- switch(xtype,
                    ppp = 1,
                    im = 1,
                    imlist = length(x))

  win <- if(nimages == 1) as.owin(x) else as.owin(x[[1]])
  win <- rescue.rectangle(win)
  rec <- as.rectangle(win)
  across <- min(diff(rec$xrange), diff(rec$yrange))
  # determine whether to expand window
  if(!expand || (6 * sig < across)) {
    result <- second.moment.engine(x, sigma=sigma, edge=edge,
                                   what=what, debug=debug, ..., varcov=varcov)
    return(result)
  }
  # need to expand window
  bigger <- grow.rectangle(rec, (7 * sig - across)/2)
  switch(xtype,
         ppp = {
           # pixellate first (to preserve pixel resolution)
           X <- pixellate(x, ..., padzero=TRUE)
           np <- npoints(x)
         },
         im = {
           X <- x
           np <- NULL
         },
         imlist = {
           X <- x
           np <- NULL
         })

  # Now expand
  if(nimages == 1) {
    X <- rebound.im(X, bigger)
    X <- na.handle.im(X, 0)
  } else {
    X <- lapply(X, rebound.im, rect=bigger)
    X <- lapply(X, na.handle.im, na.replace=0)
  }
  # Compute!
  out <- second.moment.engine(X, sigma=sigma, edge=edge,
                              what=what, debug=debug, ...,
                              obswin=win, varcov=varcov, npts=np)
  # Now clip it
  fbox <- shift(rec, origin="midpoint")
  if(nimages == 1) {
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
  } else {
    result <-
      switch(what,
             kernel     = out[fbox], 
             smooth     = lapply(out, "[", i=win),
             Kmeasure   = lapply(out, "[", i=fbox),
             Bartlett   = lapply(out, "[", i=fbox),
             edge       = out[win],
             smoothedge = list(
               smooth = lapply(out$smooth, "[", i=win),
               edge   = out$edge[win]),
             all        = list(
               kernel=out$kernel[fbox],
               smooth=lapply(out$smooth, "[", i=win),
               Kmeasure=lapply(out$Kmeasure, "[", i=fbox),
               Bartlett=lapply(out$Bartlett, "[", i=fbox),
               edge=out$edge[win]))
  }
  return(result)
}

second.moment.engine <-
  function(x, sigma=NULL, edge=TRUE,
           what=c("Kmeasure", "kernel", "smooth", 
             "Bartlett", "edge", "smoothedge", "all"),
           ...,
           kernel="gaussian",
           obswin = as.owin(x), varcov=NULL,
           npts=NULL, debug=FALSE)
{
  what <- match.arg(what)

  validate2Dkernel(kernel)
  
  if(is.ppp(x)) {
    # convert list of points to mass distribution
    X <- pixellate(x, ..., padzero=TRUE)
    if(is.null(npts))
      npts <- npoints(x)
  } else X <- x
  if(is.im(X)) {
    Xlist <- list(X)
    nimages <- 1
  } else if(all(unlist(lapply(X, is.im)))) {
    Xlist <- X
    X <- Xlist[[1]]
    nimages <- length(Xlist)
    blanklist <- vector(mode="list", length=nimages)
    names(blanklist) <- names(Xlist)
  } else stop("internal error: unrecognised format for x")
  unitsX <- unitname(X)
  # ensure obswin has same bounding frame as X
  if(!missing(obswin))
    obswin <- rebound.owin(obswin, as.rectangle(X))
  # go to work
  Y <- X$v
  Ylist <- lapply(Xlist, getElement, name="v")
  #
#  xw <- X$xrange
#  yw <- X$yrange
  # pad with zeroes
  nr <- nrow(Y)
  nc <- ncol(Y)
  Ypad <- matrix(0, ncol=2*nc, nrow=2*nr)
  Ypadlist <- rep(list(Ypad), nimages)
  for(i in 1:nimages)
    Ypadlist[[i]][1:nr, 1:nc] <- Ylist[[i]]
  Ypad <- Ypadlist[[1]]
  lengthYpad <- 4 * nc * nr
  # corresponding coordinates
#  xw.pad <- xw[1] + 2 * c(0, diff(xw))
#  yw.pad <- yw[1] + 2 * c(0, diff(yw))
  xcol.pad <- X$xcol[1] + X$xstep * (0:(2*nc-1))
  yrow.pad <- X$yrow[1] + X$ystep * (0:(2*nr-1))
  # set up kernel
  xcol.ker <- X$xstep * c(0:(nc-1),-(nc:1))
  yrow.ker <- X$ystep * c(0:(nr-1),-(nr:1))
  kerpixarea <- X$xstep * X$ystep
  if(identical(kernel, "gaussian")) {
    if(!is.null(sigma)) {
      densX.ker <- dnorm(xcol.ker, sd=sigma)
      densY.ker <- dnorm(yrow.ker, sd=sigma)
      Kern <- outer(densY.ker, densX.ker, "*") * kerpixarea
    } else if(!is.null(varcov)) {
      ## anisotropic kernel
      detSigma <- det(varcov)
      Sinv <- solve(varcov)
      halfSinv <- Sinv/2
      constker <- kerpixarea/(2 * pi * sqrt(detSigma))
      xsq <- matrix((xcol.ker^2)[col(Ypad)], ncol=2*nc, nrow=2*nr)
      ysq <- matrix((yrow.ker^2)[row(Ypad)], ncol=2*nc, nrow=2*nr)
      xy <- outer(yrow.ker, xcol.ker, "*")
      Kern <- constker * exp(-(xsq * halfSinv[1,1]
                               + xy * (halfSinv[1,2]+halfSinv[2,1])
                               + ysq * halfSinv[2,2]))
    } else 
      stop("Must specify either sigma or varcov")
  } else {
    ## non-Gaussian kernel
    ## evaluate kernel at array of points
    xker <- as.vector(xcol.ker[col(Ypad)])
    yker <- as.vector(yrow.ker[row(Ypad)])
    Kern <- evaluate2Dkernel(kernel, xker, yker,
                             sigma=sigma, varcov=varcov, ...) * kerpixarea
    Kern <- matrix(Kern, ncol=2*nc, nrow=2*nr)
  }
  # these options call for several image outputs
  if(what %in% c("all", "smoothedge"))
    result <- list()
  
  if(what %in% c("kernel", "all")) {
    # return the kernel
    # first rearrange it into spatially sensible order (monotone x and y)
    rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
    ctwist <- (-nc):(nc-1) %% (2*nc) + 1
    if(debug) {
      if(any(fave.order(xcol.ker) != rtwist))
        cat("something round the twist\n")
    }
    Kermit <- Kern[ rtwist, ctwist]
    ker <- im(Kermit, xcol.ker[ctwist], yrow.ker[ rtwist], unitname=unitsX)
    if(what == "kernel")
      return(ker)
    else 
      result$kernel <- ker
  }
  # convolve using fft
  fK <- fft(Kern)
  if(what != "edge") {
    if(nimages == 1) {
      fY <- fft(Ypad)
      sm <- fft(fY * fK, inverse=TRUE)/lengthYpad
      if(debug) {
        cat(paste("smooth: maximum imaginary part=",
                  signif(max(Im(sm)),3), "\n"))
        if(!is.null(npts))
          cat(paste("smooth: mass error=",
                    signif(sum(Mod(sm))-npts,3), "\n"))
      }
    } else {
      fYlist <- smlist <- blanklist
      for(i in 1:nimages) {
        fYlist[[i]] <- fY.i <- fft(Ypadlist[[i]])
        smlist[[i]] <- sm.i <- fft(fY.i * fK, inverse=TRUE)/lengthYpad
        if(debug) {
          cat(paste("smooth component", i, ": maximum imaginary part=",
                    signif(max(Im(sm.i)),3), "\n"))
          if(!is.null(npts))
            cat(paste("smooth component", i, ": mass error=",
                      signif(sum(Mod(sm.i))-npts,3), "\n"))
        }
      }
    }
  }
  if(what %in% c("smooth", "all", "smoothedge")) {
    # return the smoothed point pattern without edge correction
    if(nimages == 1) {
      smo <- im(Re(sm)[1:nr, 1:nc],
                xcol.pad[1:nc], yrow.pad[1:nr],
                unitname=unitsX)
      if(what == "smooth") {
        return(smo)
      } else {
        result$smooth <- smo
      }
    } else {
      smolist <- blanklist
      for(i in 1:nimages) 
        smolist[[i]] <- im(Re(smlist[[i]])[1:nr, 1:nc],
                           xcol.pad[1:nc], yrow.pad[1:nr],
                           unitname=unitsX)
      smolist <- as.solist(smolist)
      if(what == "smooth") {
        return(smolist)
      } else {
        result$smooth <- smolist
      }
    }
  }

  if(what != "edge") {
    # compute Bartlett spectrum
    if(nimages == 1) {
      bart <- BartCalc(fY, fK)  ##  bart <- Mod(fY)^2 * fK
    } else {
      bartlist <- lapply(fYlist, BartCalc, fK=fK)
    }
  }
  
  if(what %in% c("Bartlett", "all")) {
     # return Bartlett spectrum
     # rearrange into spatially sensible order (monotone x and y)
    rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
    ctwist <- (-nc):(nc-1) %% (2*nc) + 1
    if(nimages == 1) {
      Bart <- bart[ rtwist, ctwist]
      Bartlett <- im(Mod(Bart),(-nc):(nc-1), (-nr):(nr-1))
      if(what == "Bartlett")
        return(Bartlett)
      else
        result$Bartlett <- Bartlett
    } else {
      Bartlist <- blanklist
      for(i in 1:nimages) {
        Bart <- (bartlist[[i]])[ rtwist, ctwist]
        Bartlist[[i]] <- im(Mod(Bart),(-nc):(nc-1), (-nr):(nr-1))
      }
      Bartlist <- as.solist(Bartlist)
      if(what == "Bartlett")
        return(Bartlist)
      else
        result$Bartlett <- Bartlist
    }
  }
  
  #### ------- Second moment measure --------------
  #
  if(what != "edge") {
    if(nimages == 1) {
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
    } else {
      momlist <- blanklist
      for(i in 1:nimages) {
        mom.i <- fft(bartlist[[i]], inverse=TRUE)/lengthYpad
        if(debug) {
          cat(paste("2nd moment measure: maximum imaginary part=",
                    signif(max(Im(mom.i)),3), "\n"))
          if(!is.null(npts))
            cat(paste("2nd moment measure: mass error=",
                      signif(sum(Mod(mom.i))-npts^2, 3), "\n"))
        }
        mom.i <- Mod(mom.i)
        # subtract (delta_0 * kernel) * npts
        if(is.null(npts))
          stop("Internal error: second moment measure requires npts")
        mom.i <- mom.i - npts* Kern
        momlist[[i]] <- mom.i
      }
    }
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
      # apply edge correction      
      co <- fft(Mod(fM)^2 * fK, inverse=TRUE)/lengthMpad
      co <- Mod(co) 
      a <- sum(M)
      wt <- a/co
      me <- spatstat.options("maxedgewt")
      weight <- matrix(pmin.int(me, wt), ncol=2*nc, nrow=2*nr)
      # apply edge correction to second moment measure
      if(nimages == 1) {
        mom <- mom * weight
        # set to NA outside 'reasonable' region
        mom[wt > 10] <- NA
      } else {
        wgt10 <- (wt > 10)
        for(i in 1:nimages) {
          mom.i <- momlist[[i]]
          mom.i <- mom.i * weight
          # set to NA outside 'reasonable' region
          mom.i[wgt10] <- NA
          momlist[[i]] <- mom.i
        }
      }
    }
  }
  if(what != "edge") {
    # rearrange into spatially sensible order (monotone x and y)
    rtwist <- ((-nr):(nr-1)) %% (2 * nr) + 1
    ctwist <- (-nc):(nc-1) %% (2*nc) + 1
    if(nimages == 1) {
      mom <- mom[ rtwist, ctwist]
    } else {
      momlist <- lapply(momlist, "[", i=rtwist, j=ctwist)
    }
    if(debug) {
      if(any(fave.order(xcol.ker) != rtwist))
        cat("internal error: something round the twist\n")
    }
  }
  if(what %in% c("edge", "all", "smoothedge")) {
    # return convolution of window with kernel
    # (evaluated inside window only)
    con <- fft(fM * fK, inverse=TRUE)/lengthMpad
    edg <- Mod(con[1:nr, 1:nc])
    edg <- im(edg, xcol.pad[1:nc], yrow.pad[1:nr], unitname=unitsX)
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
  if(nimages == 1) {
    mom <- mom * area(obswin) / (pixarea * npts^2)
    # this is the second moment measure
    mm <- im(mom, xcol.ker[ctwist], yrow.ker[rtwist], unitname=unitsX)
    if(what == "Kmeasure")
      return(mm)
    else 
      result$Kmeasure <- mm
  } else {
    ccc <- area(obswin) / (pixarea * npts^2)
    mmlist <- blanklist
    for(i in 1:nimages) {
      mom.i <- momlist[[i]]
      mom.i <- mom.i * ccc
      # this is the second moment measure
      mmlist[[i]] <-
        im(mom.i, xcol.ker[ctwist], yrow.ker[rtwist], unitname=unitsX)
    }
    mmlist <- as.solist(mmlist)
    if(what == "Kmeasure")
      return(mmlist)
    else 
      result$Kmeasure <- mmlist
  }
  # what = "all", so return all computed objects
  return(result)
}

BartCalc <- function(fY, fK) { Mod(fY)^2 * fK }
  
Kest.fft <- function(X, sigma, r=NULL, ..., breaks=NULL) {
  verifyclass(X, "ppp")
  W <- Window(X)
  lambda <- npoints(X)/area(W)
  rmaxdefault <- rmax.rule("K", W, lambda)        
  bk <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  breaks <- bk$val
  rvalues <- bk$r
  u <- Kmeasure(X, sigma, ...)
  xx <- rasterx.im(u)
  yy <- rastery.im(u)
  rr <- sqrt(xx^2 + yy^2)
  tr <- whist(rr, breaks, u$v)
  K  <- cumsum(tr) * with(u, xstep * ystep)
  rmax <- min(rr[is.na(u$v)])
  K[rvalues >= rmax] <- NA
  result <- data.frame(r=rvalues, theo=pi * rvalues^2, border=K)
  w <- X$window
  alim <- c(0, min(diff(w$xrange), diff(w$yrange))/4)
  out <- fv(result,
            "r", quote(K(r)),
            "border",
             . ~ r, alim,
            c("r", "%s[pois](r)", "hat(%s)[fb](r)"),
            c("distance argument r",
              "theoretical Poisson %s",
              "border-corrected FFT estimate of %s"),
            fname="K",
            unitname=unitname(X)
            )
  return(out)
}

