#
#           Kmeasure.R
#
#           $Revision: 1.65 $    $Date: 2017/11/19 04:32:28 $
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
                               what=c("Kmeasure", "kernel", "smooth", 
                                 "Bartlett", "edge", "smoothedge", "all"),
                               ..., 
                               varcov=NULL, expand=FALSE, debug=FALSE) {
  if(is.null(sigma) && is.null(varcov))
    stop("must specify sigma or varcov")
  what <- match.arg(what)
  sig <- if(!is.null(sigma)) sigma else max(c(diag(varcov), sqrt(det(varcov))))

  xtype <- if(is.ppp(x)) "ppp" else
           if(is.im(x)) "im" else
           if(inherits(x, "imlist")) "imlist" else 
           if(all(sapply(x, is.im))) "imlist" else
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
  #' need to expand window
  wid <- (7 * sig - across)/2
  bigger <- grow.rectangle(rec, wid)
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

  is.second.order <- what %in% c("Kmeasure", "Bartlett", "all")
  needs.kernel <- what %in% c("kernel", "all", "Kmeasure")
  returns.several <- what %in% c("all", "smoothedge")

  # check whether Fastest Fourier Transform in the West is available
  west <- fftwAvailable()
  
  if(returns.several)
    result <- list() # several results will be returned in a list

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
  xstep <- X$xstep
  ystep <- X$ystep
  # ensure obswin has same bounding frame as X
  if(!missing(obswin))
    obswin <- rebound.owin(obswin, as.rectangle(X))
  # go to work
  Y <- X$v
  Ylist <- lapply(Xlist, getElement, name="v")
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
  xcol.pad <- X$xcol[1] + xstep * (0:(2*nc-1))
  yrow.pad <- X$yrow[1] + ystep * (0:(2*nr-1))
  # compute kernel and its Fourier transform
  if(!needs.kernel && 
     identical(kernel, "gaussian") &&
     is.numeric(sigma) && (length(sigma) == 1) &&
     spatstat.options('developer')) {
    # compute Fourier transform of kernel directly (*experimental*)
    ii <- c(0:(nr-1), nr:1)
    jj <- c(0:(nc-1), nc:1)
    zz <- -sigma^2 * pi^2/2
    uu <- exp(zz * ii^2)
    vv <- exp(zz * jj^2)
    fK <- outer(uu, vv, "*")
  } else {
    # set up kernel
    xcol.ker <- xstep * c(0:(nc-1),-(nc:1))
    yrow.ker <- ystep * c(0:(nr-1),-(nr:1))
    kerpixarea <- xstep * ystep
    if(identical(kernel, "gaussian")) {
      if(!is.null(sigma)) {
        densX.ker <- dnorm(xcol.ker, sd=sigma)
        densY.ker <- dnorm(yrow.ker, sd=sigma)
        #' WAS:  Kern <- outer(densY.ker, densX.ker, "*") * kerpixarea
        Kern <- outer(densY.ker, densX.ker, "*")
        Kern <- Kern/sum(Kern)
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
        Kern <- Kern/sum(Kern)
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
      Kern <- Kern/sum(Kern)
    }

    if(what %in% c("kernel", "all")) {
      ## kernel will be returned
      ## first rearrange it into spatially sensible order (monotone x and y)
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
    ## convolve using fft
    fK <- fft2D(Kern, west=west)
  }
  
  if(what != "edge") {
    if(nimages == 1) {
      fY <- fft2D(Ypad, west=west)
      sm <- fft2D(fY * fK, inverse=TRUE, west=west)/lengthYpad
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
        fYlist[[i]] <- fY.i <- fft2D(Ypadlist[[i]], west=west)
        smlist[[i]] <- sm.i <-
          fft2D(fY.i * fK, inverse=TRUE, west=west)/lengthYpad
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
    # compute smoothed point pattern without edge correction
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

  if(is.second.order) {
    # compute Bartlett spectrum
    if(nimages == 1) {
      bart <- BartCalc(fY, fK)  ##  bart <- Mod(fY)^2 * fK
    } else {
      bartlist <- lapply(fYlist, BartCalc, fK=fK)
    }
  }
  
  if(what %in% c("Bartlett", "all")) {
     # Bartlett spectrum will be returned
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
  if(is.second.order) {
    if(nimages == 1) {
      mom <- fft2D(bart, inverse=TRUE, west=west)/lengthYpad
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
        mom.i <- fft2D(bartlist[[i]], inverse=TRUE, west=west)/lengthYpad
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
    M <- as.mask(obswin, xy=list(x=X$xcol, y=X$yrow))$m
    # previous line ensures M has same dimensions and scale as Y 
    Mpad <- matrix(0, ncol=2*nc, nrow=2*nr)
    Mpad[1:nr, 1:nc] <- M
    lengthMpad <- 4 * nc * nr
    fM <- fft2D(Mpad, west=west)
    if(edge && is.second.order) {
      # compute kernel-smoothed set covariance
      # apply edge correction      
      co <- fft2D(Mod(fM)^2 * fK, inverse=TRUE, west=west)/lengthMpad
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
  if(is.second.order) {
    # rearrange second moment measure
    # into spatially sensible order (monotone x and y)
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
    con <- fft2D(fM * fK, inverse=TRUE, west=west)/lengthMpad
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
  pixarea <- xstep * ystep
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

