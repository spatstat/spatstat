#
#
#   Residual plots:
#         resid4plot       four panels with matching coordinates
#         resid1plot       one or more unrelated individual plots 
#         resid1panel      one panel of resid1plot
#
#   $Revision: 1.17 $    $Date: 2013/07/17 02:37:41 $
#
#

resid4plot <- function(RES, plot.neg="image", plot.smooth="imagecontour",
                       spacing=0.1, srange=NULL, monochrome=FALSE, main=NULL,
                       ...)
{
  clip     <- RES$clip
  Yclip    <- RES$Yclip
  Z        <- RES$smooth$Z
  W        <- RES$W
  Wclip    <- Yclip$window
  type     <- RES$type
  typename <- RES$typename
  Ydens    <- RES$Ydens[Wclip, drop=FALSE]
  Ymass    <- RES$Ymass[Wclip]
  # set up 2 x 2 plot with space
  wide <- diff(W$xrange)
  high <- diff(W$yrange)
  space <- spacing * max(wide,high)
  width <- wide + space + wide
  height <- high + space + high
  outerspace <- 3 * space
  plot(c(0, width) + outerspace * c(-1,1),
       c(0, height) + outerspace * c(-1,1),
       type="n", asp=1.0, axes=FALSE, xlab="", ylab="")
  # determine colour map
  if(is.null(srange)) {
    Yrange <- if(!is.null(Ydens)) summary(Ydens)$range else NULL
    Zrange <- if(!is.null(Z)) summary(Z)$range else NULL
    srange <- range(c(0, Yrange, Zrange), na.rm=TRUE)
  } else {
    stopifnot(is.numeric(srange) && length(srange) == 2)
    stopifnot(all(is.finite(srange)))
  }
    
  cols <- beachcolours(srange, if(type=="eem") 1 else 0, monochrome)
                      
  # ------ plot residuals/marks (in top left panel) ------------
  Xlowleft <- c(W$xrange[1],W$yrange[1])
  vec <- c(0, high) + c(0, space) - Xlowleft
  # shift the original window
  Ws <- shift(W, vec)
  # shift the residuals 
  Ys <- shift(Yclip,vec)
  
  # determine whether pre-plotting the window(s) is redundant
  redundant <- 
    (plot.neg == "image") && (type != "eem") && (Yclip$window$type == "mask")

  # pre-plot the window(s)
  if(!redundant) {
    if(!clip) 
      plot(Ys$window, add=TRUE, ...)
    else
      ploterodewin(Ws, Ys$window, add=TRUE, ...)
  }

  switch(plot.neg,
         discrete={
           neg <- (Ys$marks < 0)
           if(any(c("maxsize", "markscale") %in% names(list(...))))
             plot(Ys[neg], add=TRUE, ...)
           else {
             hackmax <- 0.5 * sqrt(area.owin(Wclip)/Yclip$n)
             plot(Ys[neg], add=TRUE, maxsize=hackmax, ...)
           }
           plot(Ys[!neg], add=TRUE, ...)
         },
         image={
           Yds <- shift(Ydens, vec)
           Yms <- shift(Ymass, vec)
           if(redundant)
             ploterodeimage(Ws, Yds, rangeZ=srange, colsZ=cols, ...)
           else if(type != "eem") 
             image(Yds, add=TRUE, ribbon=FALSE, col=cols, zlim=srange, ...)
           plot(Yms, add=TRUE, ...)
         }
         )
  # --------- plot smoothed surface (in bottom right panel) ------------
  vec <- c(wide, 0) + c(space, 0) - Xlowleft
  Zs <- shift.im(Z, vec)
  switch(plot.smooth,
         image={
           image(Zs, add=TRUE, col=cols, zlim=srange, ribbon=FALSE, ...)},
         contour={contour(Zs, add=TRUE, ...)},
         persp={ warning("persp not available in 4-panel plot") },
         imagecontour={
             image(Zs, add=TRUE, col=cols, zlim=srange, ribbon=FALSE, ...)
             contour(Zs, add=TRUE, ...)
           }
         )  
  lines(Zs$xrange[c(1,2,2,1,1)], Zs$yrange[c(1,1,2,2,1)])
  # -------------- lurking variable plots -----------------------
  do.lines <-
    function(x, y, defaulty=1, ...) {
      do.call("lines",
              resolve.defaults(list(x, y),
                               list(...),
                               list(lty=defaulty)))
    }
  # --------- lurking variable plot for x coordinate ------------------
  #           (cumulative or marginal)
  #           in bottom left panel
  if(!is.null(RES$xmargin)) {
    a <- RES$xmargin
    observedV <-    a$xZ
    observedX <-    a$x
    theoreticalV <- a$ExZ
    theoreticalX <- a$x
    theoreticalSD <- NULL
    ylabel <- paste("marginal of", typename)
  } else if(!is.null(RES$xcumul)) {
    a <- RES$xcumul
    observedX <- a$empirical$covariate
    observedV <- a$empirical$value
    theoreticalX <- a$theoretical$covariate
    theoreticalV <- a$theoretical$mean
    theoreticalSD <- a$theoretical$sd
    ylabel <- paste("cumulative sum of", typename)
  }
  # pretty axis marks
  pX <- pretty(theoreticalX)
  if(is.null(theoreticalSD))
    pV <- pretty(c(0,observedV,theoreticalV))
  else
    pV <- pretty(c(0,observedV,theoreticalV,
                   theoreticalV+2*theoreticalSD,
                   theoreticalV-2*theoreticalSD))
  # rescale smoothed values
  rr <- range(c(0, observedV, theoreticalV, pV))
  yscale <- function(y) { high * (y - rr[1])/diff(rr) }
  xscale <- function(x) { x - W$xrange[1] }
  do.lines(xscale(observedX), yscale(observedV), 1, ...)
  do.lines(xscale(theoreticalX), yscale(theoreticalV), 2, ...)
  if(!is.null(theoreticalSD)) {
    do.lines(xscale(theoreticalX),
             yscale(theoreticalV + 2 * theoreticalSD),
             3, ...)
    do.lines(xscale(theoreticalX),
             yscale(theoreticalV - 2 * theoreticalSD),
             3, ...)
  }
  axis(side=1, pos=0, at=xscale(pX), labels=pX)
  text(xscale(mean(theoreticalX)), - outerspace, "x coordinate")
  axis(side=2, pos=0, at=yscale(pV), labels=pV)
  text(-outerspace, yscale(mean(pV)), ylabel, srt=90)
  
  # --------- lurking variable plot for y coordinate ------------------
  #           (cumulative or marginal)
  #           in top right panel
  if(!is.null(RES$ymargin)) {
    a <- RES$ymargin
    observedV <-    a$yZ
    observedY <-    a$y
    theoreticalV <- a$EyZ
    theoreticalY <- a$y
    theoreticalSD <- NULL
    ylabel <- paste("marginal of", typename)
  } else if(!is.null(RES$ycumul)) {
    a <- RES$ycumul
    observedV <- a$empirical$value
    observedY <- a$empirical$covariate
    theoreticalY <- a$theoretical$covariate
    theoreticalV <- a$theoretical$mean
    theoreticalSD <- a$theoretical$sd
    ylabel <- paste("cumulative sum of", typename)
  }
  # pretty axis marks
  pY <- pretty(theoreticalY)
  if(is.null(theoreticalSD))
    pV <- pretty(c(0,observedV,theoreticalV))
  else
    pV <- pretty(c(0,observedV,theoreticalV,
                   theoreticalV+2*theoreticalSD,
                   theoreticalV-2*theoreticalSD))
  # rescale smoothed values
  rr <- range(c(0, observedV, theoreticalV, pV))
  yscale <- function(y) { y - W$yrange[1] + high + space}
  xscale <- function(x) { wide + space + wide * (rr[2] - x)/diff(rr) }
  do.lines(xscale(observedV), yscale(observedY), 1, ...)
  do.lines(xscale(theoreticalV), yscale(theoreticalY), 2, ...)
  if(!is.null(theoreticalSD)) {
    do.lines(xscale(theoreticalV+2*theoreticalSD),
             yscale(theoreticalY),
             3, ...)
    do.lines(xscale(theoreticalV-2*theoreticalSD),
             yscale(theoreticalY),
             3, ...)
  }
  axis(side=4, pos=width, at=yscale(pY), labels=pY)
  text(width + outerspace, yscale(mean(theoreticalY)), "y coordinate", srt=90)
  axis(side=3, pos=height, at=xscale(pV), labels=pV)
  text(xscale(mean(pV)), height + outerspace, ylabel)
  #
  if(!is.null(main))
    title(main=main)
  invisible(NULL)
}

#
#
#   Residual plot: single panel(s)
#
#

resid1plot <- function(RES, opt,
                       plot.neg="image", plot.smooth="imagecontour",
                       srange=NULL, monochrome=FALSE, main=NULL,
                       ...) {
  clip  <- RES$clip
  Y     <- RES$Y
  Yclip <- RES$Yclip
  Z     <- RES$smooth$Z
  W     <- RES$W
  Wclip <- Yclip$window
  type  <- RES$type
  Ydens <- RES$Ydens[Wclip, drop=FALSE]
  Ymass <- RES$Ymass[Wclip]
  # determine colour map
  if(opt$all || opt$marks || opt$smooth) {
    if(is.null(srange)) {
      Yrange <- if(!is.null(Ydens)) summary(Ydens)$range else NULL
      Zrange <- if(!is.null(Z)) summary(Z)$range else NULL
      srange <- range(c(0, Yrange, Zrange), na.rm=TRUE)
    }
    cols <- beachcolours(srange, if(type=="eem") 1 else 0, monochrome)
  }
  # determine main heading
  if(is.null(main)) {
    prefix <- if(opt$marks) NULL 
              else if(opt$smooth) "Smoothed"
              else if(opt$xcumul) "Lurking variable plot for x coordinate\n"
              else if(opt$ycumul) "Lurking variable plot for y coordinate\n"
              else if(opt$xmargin) "Lurking variable plot for x coordinate\n"
              else if(opt$ymargin) "Lurking variable plot for y coordinate\n"
    main <- paste(prefix, RES$typename)
  }
    
  # ------------- residuals ---------------------------------
  if(opt$marks) {
    # determine whether pre-plotting the window(s) is redundant
    redundant <- 
      (plot.neg == "image") && (type != "eem") && (Yclip$window$type == "mask")
  # pre-plot the window(s)
    if(redundant)
      plot(as.rectangle(W), box=FALSE, main="", ...)
    else {
      if(!clip) 
        plot(W, main="", ...)
      else
        ploterodewin(W, Wclip, main="", ...)
    }

    switch(plot.neg,
           discrete={
             neg <- (Y$marks < 0)
             if(any(c("maxsize", "markscale") %in% names(list(...))))
               plot(Y[neg], add=TRUE, ...)
             else {
               hackmax <- 0.5 * sqrt(area.owin(Wclip)/Yclip$n)
               plot(Y[neg], add=TRUE, maxsize=hackmax, ...)
             }
             plot(Y[!neg], add=TRUE, ...)
           },
         image={
           if(redundant)
             ploterodeimage(W, Ydens, rangeZ=srange, colsZ=cols, ...)
           else if(type != "eem") 
             image(Ydens, col=cols, zlim=srange, add=TRUE, ribbon=FALSE, ...)
           plot(Ymass, add=TRUE, ...)
         }
         )
    title(main=main)
  }
  # -------------  smooth -------------------------------------
  if(opt$smooth) {
    if(!clip) {
      switch(plot.smooth,
           image={image(Z, main=main, axes=FALSE, xlab="", ylab="",
                        col=cols, zlim=srange, ribbon=FALSE, ...)},
           contour={contour(Z, main=main, axes=FALSE, xlab="", ylab="", ...)},
           persp={persp(Z, main=main, axes=FALSE, xlab="", ylab="", ...)},
           imagecontour={
             image(Z, main=main, axes=FALSE, xlab="", ylab="",
                   col=cols, zlim=srange, ribbon=FALSE, ...)
             contour(Z, add=TRUE, ...)
           }
           )
    } else {
      switch(plot.smooth,
           image={
             plot(as.rectangle(W), box=FALSE, main=main, ...)
             ploterodeimage(W, Z, colsZ=cols, rangeZ=srange, ...)
           },
           contour={
             plot(W, main=main, ...)
             contour(Z, add=TRUE, ...)
           },
           persp={
             persp(Z, main=main, axes=FALSE, xlab="", ylab="", ...)
           # there is no 'add' option for 'persp'
           },
           imagecontour={
             plot(as.rectangle(W), box=FALSE, main=main, ...)
             ploterodeimage(W, Z, colsZ=cols, rangeZ=srange, ...)
             contour(Z, add=TRUE, ...)
           }
           )
    }
  }
  # ------------  cumulative x -----------------------------------------
  if(opt$xcumul) {
    a <- RES$xcumul
    obs <- a$empirical
    theo <- a$theoretical
    resid1panel(obs$covariate, obs$value,
               theo$covariate, theo$mean, theo$sd,
               "x coordinate", "cumulative mark", main=main, ...)
  }
  
  # ------------  cumulative y -----------------------------------------
  if(opt$ycumul) {
    a <- RES$ycumul
    obs <- a$empirical
    theo <- a$theoretical
    resid1panel(obs$covariate, obs$value,
               theo$covariate, theo$mean, theo$sd,
               "y coordinate", "cumulative mark", main=main, ...)
  }
  # ------------  x margin -----------------------------------------
  if(opt$xmargin) {
    a <- RES$xmargin
    resid1panel(a$x, a$xZ, a$x, a$ExZ, NULL,
               "x coordinate", "marginal of residuals", main=main, ...)
  }
  # ------------  y margin -----------------------------------------
  if(opt$ymargin) {
    a <- RES$ymargin
    resid1panel(a$y, a$yZ, a$y, a$EyZ, NULL,
               "y coordinate", "marginal of residuals", main=main, ...)
  }
  
  return(invisible(NULL))
}


resid1panel <- function(observedX, observedV,
                        theoreticalX, theoreticalV, theoreticalSD, xlab, ylab,
                        ...)
{
    # work out plot range
    rX <- range(observedX, theoreticalX)
    rV <- range(c(0, observedV, theoreticalV))
    if(!is.null(theoreticalSD))
        rV <- range(c(rV, theoreticalV + 2*theoreticalSD,
                          theoreticalV - 2*theoreticalSD))
    # argument handling
    do.lines <-
      function(x, y, defaulty=1, ...) {
        do.call("lines",
                resolve.defaults(list(x, y),
                                 list(...),
                                 list(lty=defaulty)))
      }
    # start plot
    plot(rX, rV, type="n", xlab=xlab, ylab=ylab, ...)
    do.lines(observedX, observedV, 1, ...)
    do.lines(theoreticalX, theoreticalV, 2, ...)
    if(!is.null(theoreticalSD)) {
      do.lines(theoreticalX, theoreticalV + 2 * theoreticalSD, 3, ...)
      do.lines(theoreticalX, theoreticalV - 2 * theoreticalSD, 3, ...)
    }
}

#
#
ploterodewin <- function(W1, W2, col.edge=grey(0.75), col.inside=rgb(1,0,0),
                         ...) {
  # internal use only
  # W2 is assumed to be an erosion of W1
  switch(W1$type,
         rectangle={
           plot(W1, ...)
           plot(W2, add=TRUE, lty=2)
         },
         polygonal={
           plot(W1, ...)
           plot(W2, add=TRUE, lty=2)
         },
         mask={
           Z <- as.im(W1)
           x <- as.vector(raster.x(W1))
           y <- as.vector(raster.y(W1))
           ok <- inside.owin(x, y, W2)
           Z$v[ok] <- 2
           plot(Z, ..., col=c(col.edge, col.inside),
                  add=TRUE, ribbon=FALSE)
         }
         )
}

ploterodeimage <- function(W, Z, ..., Wcol=grey(0.75), rangeZ, colsZ) {
  # Internal use only
  # Image Z is assumed to live on a subset of mask W
  # colsZ are the colours for the values in the range 'rangeZ'

  if(W$type != "mask") {
    plot(W, add=TRUE)
    W <- as.mask(W)
  }
  
  # Extend the colour map to include an extra colour for pixels in W
  # (1) Add the desired colour of W to the colour map
  pseudocols <- c(Wcol, colsZ)
  # (2) Breakpoints
  bks <- seq(from=rangeZ[1], to=rangeZ[2], length=length(colsZ)+1)
  dZ <- diff(bks)[1]
  pseudobreaks <- c(rangeZ[1] - dZ, bks)
  # (3) Determine a fake value for pixels in W
  Wvalue <- rangeZ[1] - dZ/2

  # Create composite image on W grid
  # (with W-pixels initialised to Wvalue)
  X <- as.im(Wvalue, W)
  # Look up Z-values of W-pixels
  xx <- as.vector(raster.x(W))
  yy <- as.vector(raster.y(W))
  Zvalues <- lookup.im(Z, xx, yy, naok = TRUE, strict=FALSE)
  # Overwrite pixels in Z
  inZ <- !is.na(Zvalues)
  X$v[inZ] <- Zvalues[inZ]

  image(X, ..., add=TRUE, ribbon=FALSE, col=pseudocols, breaks=pseudobreaks)
  return(list(X, pseudocols, pseudobreaks))
}


  
