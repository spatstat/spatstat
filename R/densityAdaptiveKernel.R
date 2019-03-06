#'
#'   densityAdaptiveKernel.R
#'
#'   $Revision: 1.3 $  $Date: 2019/03/06 02:54:35 $
#'
#'
#'  Adaptive kernel smoothing via 3D FFT
#'

densityAdaptiveKernel <- function(X, ...) {
  UseMethod("densityAdaptiveKernel")
}

densityAdaptiveKernel.ppp <- local({

  ada <- function(X, bw,
                  ...,
                  bwrange=range(bw), 
                  edge=TRUE, diggle=TRUE, weights=NULL, ngroups=NULL,
                  nbw=64, what=c("smooth", "all"), ncoarse=NULL,
                  taperfac=c(5,10), verbose=TRUE) {
    stopifnot(is.ppp(X))
    nX <- npoints(X)
    what <- match.arg(what)
  
    check.nvector(taperfac,2,things="taper functionals")
    if(any(taperfac<=0)) stop("'taperfac' components must be positive")

    if(missing(bw) || is.null(bw)) 
      bw <- do.call.matched(bw.abram,
                            resolve.defaults(list(X=X, at="pixels"),
                                             list(...)),
                            extrargs=names(args(as.mask)))

    if(is.numeric(bw)) {
      check.nvector(bw, nX, oneok=TRUE)
      if(length(bw) == 1) bw <- rep(bw, nX)
      bwX <- bw
    } else if(is.im(bw)) {
      bwX <- bw[X, drop=FALSE]
      if(anyNA(bwX))
        stop("Some data points lie outside the domain of image 'bw'",
             call.=FALSE)
    } else stop("Argument 'bw' should be a numeric vector or a pixel image")

    check.range(bwrange)

    marks(X) <- bwX
  
    if(edge && !diggle && !is.im(bw))
      stop("For the uniform edge correction, 'bw' must be a pixel image",
           call.=FALSE)
    
    if(weighted <- !is.null(weights)) {
      check.nvector(weights, nX, oneok=TRUE)
      if(length(weights) == 1) weights <- rep(weights, nX)
    } else weights <- rep(1,nX)
    
    splitting <- !is.null(ngroups)
    if(splitting) {
      check.1.integer(ngroups)
      stopifnot(ngroups > 0)
    }
  
    if(verbose && !splitting) cat("Discretising...")
  
    #' discretise window
    WX <- Window(X)
    isrect <- is.rectangle(WX)
    W <- as.mask(WX, ...)
    insideW <- W$m
    dimW <- W$dim
    nr <- dimW[1]
    nc <- dimW[2]
    xstep <- W$xstep
    ystep <- W$ystep
   
    if(splitting){
      nbwq <- ceiling(nbw/ngroups)
      bwquan <- quantile(bwX,seq(0,1,1/ngroups))
      bwqpts <- findInterval(bwX,bwquan,all.inside=TRUE)
	
      smosum <- 0

      for(hq in 1:ngroups){
        if(verbose) cat("Partition",hq,"of",ngroups,"\b...\n")
        Xq <- X[bwqpts==hq]
        bwq <- marks(Xq)
	     
        weightsq <- if(weighted) weights[bwqpts==hq] else NULL
        ies <- if(diggle&&edge) ncoarse else NULL
	   
        smobit <- densityAdaptiveKernel.ppp(X=Xq, bw=bwq,
                                            edge = edge && diggle,
                                            diggle=diggle,
                                            weights=weightsq,
                                            what="smooth",
                                            nbw=nbwq, ncoarse=ies,
                                            verbose=FALSE, 
                                            taperfac=taperfac, ...)
        smosum <- smosum + smobit
      }
      if(!edge || diggle)
        return(smosum)
      #' uniform edge correction is deferred until further below..
    }
  
    #' discretise spatial locations of data points
    ij <- nearest.raster.point(X$x, X$y, W)
  
    #' discretise bandwidth values
    zrange <- range(0, 1.5 * bwrange)
    zbreaks <- seq(zrange[1], zrange[2], length=nbw+1)
    zstep <- diff(zbreaks[1:2])
    zvalues <- zbreaks[-1] - zstep/2
    kslice <- findInterval(bwX, zbreaks, all.inside=TRUE)

    #' grid coordinates (padded)
    xcol.pad <- W$xcol[1] + xstep * (0:(2*nc-1))
    yrow.pad <- W$yrow[1] + ystep * (0:(2*nr-1))
    z.pad    <- zvalues[1] + zstep * (0:(2*nbw - 1))
    
    if(verbose) cat("Done.\n Forming kernel...")
  
    ## set up kernel
    xcol.ker <- xstep * c(0:(nc-1),-(nc:1))
    yrow.ker <- ystep * c(0:(nr-1),-(nr:1))
    z.ker    <- zstep   * c(0:(nbw-1), -(nbw:1))
    pixarea <- xstep * ystep
    kerpixvol <- xstep * ystep * zstep
    az <- abs(z.ker)

    mapback <- function(scope, a) { (a - scope[1])/diff(scope) }
    limL <- mapback(c(-2, taperfac[1]), c(-3,3)) * bwrange[1]
    limR <- mapback(c(-taperfac[2], 4), c(-3,3)) / bwrange[2]
    taper <- taperoff(az, limL[1], limL[2]) * taperoff(1/az, limR[1], limR[2])
    ## taper <- ifelse(az < bwrange[1],
    ##                  norm.decay(az,bwrange[1],c(-2,taperfac[1]),"left"),#'ALTERNATIVE:  exp(taperfac[1] * (az/bwrange[1] - 1)),
    ##                  ifelse(az > bwrange[2],
    ##                         norm.decay(az,bwrange[2],c(-taperfac[2],4),"right"),#'ALTERNATIVE: exp(-taperfac[2] * (az/bwrange[2] - 1)), 
    ##                         1))
    Kern <- array(0, dim=2*c(nr, nc, nbw))
    for(k in 1:(2*nbw)) {
      azk <- az[k]
      if(azk > 0) {
        densX.ker <- dnorm(xcol.ker, sd=azk)
        densY.ker <- dnorm(yrow.ker, sd=azk)
        Kern[,,k] <- outer(densY.ker, densX.ker, "*") * taper[k] * pixarea
      }
    }
    if(verbose) cat("Done.\n Taking FFT of kernel...")
    ## Fourier transform of kernel
    fK <- fft(Kern)
    if(verbose) cat("Done.\n")

    pixarea_sub <- WK <- NULL
    edgeW <- 1
    edgeX <- pixarea
    if(edge) {
      if(is.null(ncoarse)) {
        #' convolve window with kernel
        Wpad <- array(0, dim=2*c(nr, nc, nbw))
        Wpad[1:nr, 1:nc, 1] <- W$m * pixarea
        if(verbose) cat("FFT of window...")
        fW <- fft(Wpad)
        if(verbose) cat("Inverse FFT of smoothed window...")
        WK <- fft(fW * fK, inverse=TRUE)/prod(dim(Wpad))
        if(verbose)
          cat(paste("\n[Window convolution: maximum imaginary part=",
                    signif(max(abs(Im(WK))),3), "]\n"))
        if(verbose) cat("Looking up edge correction weights...")
        if(diggle) {
          #' Diggle-type edge correction:
          #' edge weights associated with data points
          #' look up edge weights
          edgeX <- Re(WK)[cbind(ij$row, ij$col, kslice)]
        } else {
          #' Uniform-type edge correction: edge weights associated with pixels
          edgeW <- bwW <- bw[W, drop=FALSE]
          kim <- eval.im(findInterval(bwW, zbreaks, all.inside=TRUE))
          iim <- as.im(row(bwW), W=bwW)
          jim <- as.im(col(bwW), W=bwW)
          df <- pairs(iim, jim, kim, plot=FALSE)
          edgeW[] <- Re(WK)[as.matrix(df)]/pixarea
        }
        weights <- weights/edgeX
        if(verbose) cat("Done\n")
      } else {
        if(!diggle){
          edgeW <- do_interp(ncoarse,WX,X,zrange,bw,bwrange,pixarea,weights,taperfac,FALSE)
          weights <- weights/edgeX
        } else {
          weights <- do_interp(ncoarse,WX,X,zrange,bwX,bwrange,pixarea,weights,taperfac,TRUE)
          edgeX <- 1/weights  
        }
      }
    } else {
      weights <- weights/edgeX
    }

    if(splitting && edge && !diggle) {
      #' convolution has already been computed: apply uniform edge correction 
      smo <- (smosum/edgeW)[W, drop=FALSE]
    } else {
      #' Convolution with data points
    
      #' generate table of 3D points (with padding)
      if(verbose) cat("Discretising point locations...")
 
      rowfac <- factor(ij$row, levels=1:(2*nr))
      colfac <- factor(ij$col, levels=1:(2*nc))
      kfac <- factor(kslice, levels=1:(2*nbw))

      Xpad <- tapplysum(weights, list(rowfac, colfac, kfac))
  
      Xpad <- unname(unclass(Xpad))
      #' convolve point masses with kernel
      if(verbose) cat("Done.\nFFT of point locations...")
      fX <- fft(Xpad)
      if(verbose) cat("Inverse FFT of smoothed point locations...")
      sm <- fft(fX * fK, inverse=TRUE)/prod(dim(Xpad))
      if(verbose){
        cat("Done.\n")
        cat(paste("maximum imaginary part=",
                  signif(max(abs(Im(sm))),3), "\n"))
      }
      #' slice kernel at z=0
      smo <- im(Re(sm)[1:nr, 1:nc, 1],
                xcol.pad[1:nc], yrow.pad[1:nr],
                unitname=unitname(X))
    }
  
    #' clip to original window
    if(!isrect)
      smo <- smo[WX, drop=FALSE]
    #' return
    if(what == "smooth")
      return(smo/edgeW)
    everything <- list(smooth=smo/edgeW,
                       Xpad=Xpad, W=W, Kern=Kern, 
                       xcol.ker=xcol.ker, yrow.ker=yrow.ker, z.ker=z.ker,
                       WK=WK, edgeX=edgeX, edgeW=edgeW, 
                       taper=taper, kerpixvol=kerpixvol, pixarea=pixarea,
                       ps=pixarea_sub)
    return(everything)
  }

  do_interp <- function(ncoarse,WX,X,zrange,bw,bwrange,pixarea,
                        weights,taperfac,diggle){
    nbw_sub <- ncoarse[3]
    W_sub <- as.mask(WX,dimyx=ncoarse[2:1])
    dimW_sub <- W_sub$dim
    nr_sub <- dimW_sub[1]
    nc_sub <- dimW_sub[2]
    xstep_sub <- W_sub$xstep
    ystep_sub <- W_sub$ystep
    ij_sub <- nearest.raster.point(X$x, X$y, W_sub)
    zbreaks_sub <- seq(zrange[1], zrange[2], length=nbw_sub+1)
    zstep_sub <- diff(zbreaks_sub[1:2])
    zvalues_sub <- zbreaks_sub[-1] - zstep_sub/2
    kslice_sub <- findInterval(bw, zbreaks_sub, all.inside=TRUE)
    xcol.ker_sub <- xstep_sub * c(0:(nc_sub-1),-(nc_sub:1))
    yrow.ker_sub <- ystep_sub * c(0:(nr_sub-1),-(nr_sub:1))
    z.ker_sub    <- zstep_sub   * c(0:(nbw_sub-1), -(nbw_sub:1))
    pixarea_sub <- xstep_sub * ystep_sub
    kerpixvol_sub <- xstep_sub * ystep_sub * zstep_sub
    az_sub <- abs(z.ker_sub)
    limL <- c(-3,3) * bwrange[1]/(taperfac[1] + 2)
    limR <- c(-3,3) * (1/bwrange[2])/(taperfac[2] + 4)
    taper <- taperoff(az_sub, limL[1], limL[2]) * taperoff(1/az_sub, limR[1], limR[2])
#  taper_sub <- ifelse(az_sub < bwrange[1],
#                      norm.decay(az_sub,bwrange[1],c(-2,taperfac[1]),"left"), #'ALTERNATIVE: exp(taperfac[1] * (az/bwrange[1] - 1)),
#                      ifelse(az_sub > bwrange[2],
#                             norm.decay(az_sub,bwrange[2],c(-taperfac[2],4),"right"), #'ALTERNATIVE: exp(-taperfac[2] * (az/bwrange[2] - 1)),
#                             1))
    Kern_sub <- array(0, dim=2*c(nr_sub, nc_sub, nbw_sub))
    for(k in 1:(2*nbw_sub)) {
      azk <- az_sub[k]
      if(azk > 0) {
        densX.ker_sub <- dnorm(xcol.ker_sub, sd=azk)
        densY.ker_sub <- dnorm(yrow.ker_sub, sd=azk)
        Kern_sub[,,k] <- outer(densY.ker_sub, densX.ker_sub, "*") * taper_sub[k] * pixarea_sub
      }
    }
    fK_sub <- fft(Kern_sub)
    Wpad_sub <- array(0, dim=2*c(nr_sub, nc_sub, nbw_sub))
    Wpad_sub[1:nr_sub, 1:nc_sub, 1] <- W_sub$m * pixarea_sub
    fW_sub <- fft(Wpad_sub)
    WK_sub <- fft(fW_sub * fK_sub, inverse=TRUE)/prod(dim(Wpad_sub))
    
    WK_grint <- Re(WK_sub)[1:ncoarse[1],1:ncoarse[2],1:ncoarse[3]]
    dimnames(WK_grint) <- list(W_sub$yrow,W_sub$xcol,zvalues_sub)
  
    if(diggle){
      X_grint <- cbind(X$y, X$x, bw)
      WK_intrp <- interpolate(X_grint,WK_grint)/pixarea_sub*pixarea
      weights <- (weights/WK_intrp)
      return(weights)
    } else {
      WK_ims <- apply(WK_grint,3,function(x) as.im(interp.im,W=WX,Z=im(x,xcol=W_sub$xcol,yrow=W_sub$yrow),dimyx=dim(bw)))
      WK_intrp <- as.matrix(WK_ims[[1]])
      for(i in 2:length(WK_ims))
        WK_intrp <- abind(WK_intrp,as.matrix(WK_ims[[i]]),along=3)
    
      edgeW_int <- bwW_int <- bw[WX, drop=FALSE]
      kim_sub <- eval.im(findInterval(bwW_int, zbreaks_sub, all.inside=TRUE))
      iim_sub <- as.im(row(bwW_int), W=bwW_int)
      jim_sub <- as.im(col(bwW_int), W=bwW_int)
      df_sub <- pairs(iim_sub, jim_sub, kim_sub, plot=FALSE)
      edgeW_int[] <- Re(WK_intrp)[as.matrix(df_sub)]/pixarea_sub
      return(edgeW_int)
    }
  }

  ada
})
