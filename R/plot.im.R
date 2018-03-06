#
#   plot.im.R
#
#  $Revision: 1.122 $   $Date: 2017/10/11 04:33:39 $
#
#  Plotting code for pixel images
#
#  plot.im
#  image.im
#  contour.im
#
###########################################################################

plot.im <- local({

  ## auxiliary functions

  image.doit <- function(imagedata, ...,
                         extrargs=graphicsPars("image"), W,
                         workaround=FALSE) {
    aarg <- resolve.defaults(...)
    add      <- resolve.1.default(list(add=FALSE),     aarg)
    show.all <- resolve.1.default(list(show.all=!add), aarg)
    addcontour <- resolve.1.default(list(addcontour=FALSE), aarg)
    args.contour <- resolve.1.default(list(args.contour=list()), aarg)
    if(add && show.all) {
      ## set up the window space *with* the main title
      ## using the same code as plot.owin, for consistency
      do.call.matched(plot.owin,
                      resolve.defaults(list(x=W, type="n"), aarg), 
                      extrargs=graphicsPars("owin"))
    }
    if(workaround && identical(aarg$useRaster, TRUE)) {
      #' workaround for bug 16035
      #' detect reversed coordinates
      usr <- par('usr')
      xrev <- (diff(usr[1:2]) < 0) 
      yrev <- (diff(usr[3:4]) < 0)
      if(xrev || yrev) {
        #' flip matrix of pixel values, because the device driver does not
        z <- imagedata$z
        d <- dim(z) # z is in the orientation expected for image.default
        if(xrev) z <- z[d[1]:1,       , drop=FALSE]
        if(yrev) z <- z[      , d[2]:1, drop=FALSE]
        imagedata$z <- z
      }
    }
    extrargs <- setdiff(extrargs, c("claim.title.space", "box"))
    z <- do.call.matched(image.default,
                         append(imagedata, aarg),
                         extrargs=extrargs)
    if(addcontour)
      do.call(do.contour,
              resolve.defaults(imagedata,
                               list(add=TRUE),
                               args.contour,
                               list(col=par('fg')),
                               aarg,
                               .StripNull=TRUE))
    return(z)
  }

  do.contour <- function(x, y, z, ..., drawlabels=TRUE) {
    nx <- length(x)
    ny <- length(y)
    nz <- dim(z)
    if(nx > nz[1]) {
      if(nz[1] == 1) {
        z <- rbind(z, z)
        nz <- dim(z)
        drawlabels <- FALSE
      } else {
        x <- (x[-1] + x[-nx])/2
        nx <- nx-1
      }
    }
    if(ny > nz[2]) {
      if(nz[2] == 1) {
        z <- cbind(z, z)
        nz <- dim(z)
        drawlabels <- FALSE
      } else {
        y <- (y[-1] + y[-ny])/2
        ny <- ny-1
      }
    }
    do.call.matched(contour.default,
                    list(x=x, y=y, z=z, ..., drawlabels=drawlabels))
  }
                 
  do.box.etc <- function(bb, add, argh)
    do.call(box.etc, append(list(bb=bb, add=add), argh))
  
  box.etc <- function(bb, ..., add=FALSE, axes=FALSE, box=!add) {
    # axes for image
    xr <- bb$xrange
    yr <- bb$yrange
    if(box)
      rect(xr[1], yr[1], xr[2], yr[2])
    if(axes) {
      px <- pretty(xr)
      py <- pretty(yr)
      do.call.plotfun(graphics::axis,
                      resolve.defaults(
                                       list(side=1, at=px), 
                                       list(...),
                                       list(pos=yr[1])),
                      extrargs=graphicsPars("axis"))
      do.call.plotfun(graphics::axis,
                      resolve.defaults(
                                       list(side=2, at=py), 
                                       list(...),
                                       list(pos=xr[1])),
                      extrargs=graphicsPars("axis"))
    }
  }
  
  
  clamp <- function(x, v, tol=0.02 * diff(v)) {
    ok <- (x >= v[1] - tol) & (x <= v[2] + tol)
    x[ok]
  }
  
  cellbreaks <- function(x, dx) {
    nx <- length(x)
    seq(x[1] - dx/2, x[nx] + dx/2, length.out=nx+1)
  }

  log10orNA <- function(x) {
    y <- rep(NA_real_, length(x))
    ok <- !is.na(x) & (x > 0)
    y[ok] <- log10(x[ok])
    return(y)
  }

  Ticks <- function(usr, log=FALSE, nint=NULL, ...) {
    #' same as axisTicks but accepts nint=NULL as if it were missing
    if(!is.null(nint)) return(axisTicks(usr=usr, log=log, nint=nint, ...))
    return(axisTicks(usr=usr, log=log, ...))
  }
  
  # main function
  PlotIm <- function(x, ...,
                     main, 
                     add=FALSE, clipwin=NULL,
                     col=NULL, valuesAreColours=NULL, log=FALSE,
                     ncolours=256, gamma=1, 
                     ribbon=show.all, show.all=!add,
                     ribside=c("right", "left", "bottom", "top"),
                     ribsep=0.15, ribwid=0.05, ribn=1024,
                     ribscale=1, ribargs=list(), colargs=list(),
                     useRaster=NULL, workaround=FALSE,
                     do.plot=TRUE) {
    if(missing(main)) main <- short.deparse(substitute(x))
    verifyclass(x, "im")
    if(x$type == "complex") {
      cl <- match.call()
      cl$x <- solist(Re=Re(x), Im=Im(x), Mod=Mod(x), Arg=Arg(x))
      cl[[1]] <- as.name('plot')
      cl$main <- main
      out <- eval(cl, parent.frame())
      return(invisible(out))
    }
    ribside <- match.arg(ribside)
    col.given <- !is.null(col)
    dotargs <- list(...)

    stopifnot(is.list(ribargs))
    user.ticks <- ribargs$at
    user.nint <- ribargs$nint
    
    if(!is.null(clipwin)) {
      x <- x[as.rectangle(clipwin)]
      if(!is.rectangle(clipwin)) x <- x[clipwin, drop=FALSE]
    }

    zlim <- dotargs$zlim

    x <- repair.image.xycoords(x)

    xtype <- x$type
    xbox <- as.rectangle(x)
    
    do.log <- identical(log, TRUE)
    if(do.log && !(x$type %in% c("real", "integer")))
      stop(paste("Log transform is undefined for an image of type",
                 sQuote(xtype)))

    # determine whether pixel values are to be treated as colours
    if(!is.null(valuesAreColours)) {
      # argument given - validate
      stopifnot(is.logical(valuesAreColours))
      if(valuesAreColours) {
        # pixel values must be factor or character
        if(!xtype %in% c("factor", "character")) {
          warning(paste("Pixel values of type", sQuote(xtype),
                        "are not interpretable as colours"))
          valuesAreColours <- FALSE
        } else if(col.given) {
          # colour info provided: contradictory
          warning(paste("Pixel values are taken to be colour values,",
                        "because valuesAreColours=TRUE;", 
                        "the colour map (argument col) is ignored"),
                  call.=FALSE)
          col <- NULL
        }
        if(do.log) 
          warning(paste("Pixel values are taken to be colour values,",
                        "because valuesAreColours=TRUE;", 
                        "the argument log=TRUE is ignored"),
                  call.=FALSE)
      }
    } else if(col.given) {
      # argument 'col' controls colours
      valuesAreColours <- FALSE
    } else if(spatstat.options("monochrome")) {
      valuesAreColours <- FALSE
    } else {
      ## default : determine whether pixel values are colours
      strings <- switch(xtype,
                        character = { as.vector(x$v) },
                        factor    = { levels(x) },
                        { NULL })
      valuesAreColours <- is.character(strings) && 
      !inherits(try(col2rgb(strings), silent=TRUE), "try-error")
      if(valuesAreColours)
        cat("Interpreting pixel values as colours\n")
    }
    # 
    if(valuesAreColours) {
      # colour-valued images are plotted using the code for factor images
      # with the colour map equal to the levels of the factor
      switch(xtype,
             factor = {
               col <- levels(x)
             },
             character = {
               x <- eval.im(factor(x))
               xtype <- "factor"
               col <- levels(x)
             },
             {
               warning(paste("Pixel values of type", sQuote(xtype),
                             "are not interpretable as colours"))
             })
      # colours not suitable for ribbon
      ribbon <- FALSE
    } 
    
    # transform pixel values to log scale?
    if(do.log) {
      rx <- range(x, finite=TRUE)
      if(all(rx > 0)) {
        x <- eval.im(log10(x))
      } else {
        if(any(rx < 0)) 
          warning(paste("Negative pixel values",
                        "omitted from logarithmic colour map;",
                        "range of values =", prange(rx)),
                  call.=FALSE)
        if(!all(rx < 0))
          warning("Zero pixel values omitted from logarithmic colour map",
                  call.=FALSE)
        x <- eval.im(log10orNA(x))
      } 
      xtype <- x$type
      Log <- log10
      Exp <- function(x) { 10^x }
    } else {
      Log <- Exp <- function(x) { x }
    }
    
    imagebreaks <- NULL
#    ribbonvalues <- ribbonbreaks <- NULL
    ribbonvalues <- NULL

    ## NOW DETERMINE THE COLOUR MAP
    colfun <- colmap <- NULL
    if(valuesAreColours) {
      ## pixel values are colours; set of colours was determined earlier
      colmap <- colourmap(col=col, inputs=col)
    } else if(!col.given) {
      ## no colour information given: use default
      colfun <- spatstat.options("image.colfun")
    } else if(inherits(col, "colourmap")) {
      ## Bob's your uncle
      colmap <- col
    } else if(is.function(col)) {
      ## Some kind of function determining a colour map
      if(names(formals(col))[1] == "n") {
        ## function(n) -> colour values
        colfun <- col
      } else {
        ## colour map determined by a rule (e.g. 'beachcolours')
        colmap <- invokeColourmapRule(col, x, zlim=zlim, colargs=colargs)
        if(is.null(colmap))
          stop("Unrecognised syntax for colour function")
      }
    }

    switch(xtype,
           real    = {
             vrange <- range(x, finite=TRUE)
             vrange <- range(zlim, vrange)
             if(!is.null(colmap)) {
               # explicit colour map
               s <- summary(colmap)
               if(s$discrete)
                 stop("Discrete colour map is not applicable to real values")
               imagebreaks <- s$breaks
               vrange <- range(imagebreaks)
               col <- s$outputs
             } 
             trivial <- (diff(vrange) <= .Machine$double.eps)
             if(!trivial) {
               # ribbonvalues: domain of colour map (pixel values)
               # ribbonrange: (min, max) of pixel values in image
               # nominalrange: range of values shown on ribbon 
               # nominalmarks: values shown on ribbon at tick marks
               # ribbonticks: pixel values of tick marks 
               # ribbonlabels: text displayed at tick marks
               ribbonvalues <- seq(from=vrange[1], to=vrange[2],
                                   length.out=ribn)
               ribbonrange <- vrange
               nominalrange <- Log(ribscale * Exp(ribbonrange))
               nominalmarks <- user.ticks %orifnull% Ticks(nominalrange,
                                                           log=do.log,
                                                           nint=user.nint)
               ribbonticks <- Log(nominalmarks/ribscale)
               ribbonlabels <- paste(nominalmarks)
             }
           },
           integer = {
             values <- as.vector(x$v)
             values <- values[!is.na(values)]
             uv <- unique(values)
             vrange <- range(uv, finite=TRUE)
             vrange <- range(zlim, vrange)
             nvalues <- length(uv)
             trivial <- (nvalues < 2)
             if(!trivial){
               nominalrange <- Log(ribscale * Exp(vrange))
               if(!is.null(user.ticks)) {
                 nominalmarks <- user.ticks
               } else {
                 nominalmarks <- Ticks(nominalrange,
                                       log=do.log,
                                       nint = user.nint)
                 nominalmarks <- nominalmarks[nominalmarks %% 1 == 0]
               }
               ribbonticks <- Log(nominalmarks/ribscale)
               ribbonlabels <- paste(nominalmarks)
               if(!do.log && identical(all.equal(ribbonticks,
                                                 vrange[1]:vrange[2]), TRUE)) {
                 # each possible pixel value will appear in ribbon
                 ribbonvalues <- vrange[1]:vrange[2]
                 imagebreaks <- c(ribbonvalues - 0.5, vrange[2] + 0.5)
                 ribbonrange <- range(imagebreaks)
                 ribbonticks <- ribbonvalues
                 ribbonlabels <- paste(ribbonticks * ribscale)
               } else {
                 # not all possible values will appear in ribbon
                 ribn <- min(ribn, diff(vrange)+1)
                 ribbonvalues <- seq(from=vrange[1], to=vrange[2],
                                     length.out=ribn)
                 ribbonrange <- vrange
               }
             }
             if(!is.null(colmap)) {
               # explicit colour map
               s <- summary(colmap)
               imagebreaks <-
                 if(!s$discrete) s$breaks else
                 c(s$inputs[1] - 0.5, s$inputs + 0.5)
               col <- s$outputs
             }
           },
           logical = {
             values <- as.integer(as.vector(x$v))
             values <- values[!is.na(values)]
             uv <- unique(values)
             trivial <- (length(uv) < 2)
             vrange <- c(0,1)
             imagebreaks <- c(-0.5, 0.5, 1.5)
             ribbonvalues <- c(0,1)
             ribbonrange <- range(imagebreaks)
#             ribbonbreaks <- imagebreaks
             ribbonticks <- user.ticks %orifnull% ribbonvalues
             ribbonlabels <- c("FALSE", "TRUE")
             if(!is.null(colmap)) 
               col <- colmap(c(FALSE,TRUE))
           },
           factor  = {
             lev <- levels(x)
             nvalues <- length(lev)
             trivial <- (nvalues < 2)
             # ensure all factor levels plotted separately
             fac <- factor(lev, levels=lev)
             intlev <- as.integer(fac)
             imagebreaks <- c(intlev - 0.5, max(intlev) + 0.5)
             ribbonvalues <- intlev
             ribbonrange <- range(imagebreaks)
#             ribbonbreaks <- imagebreaks
             ribbonticks <- user.ticks %orifnull% ribbonvalues
             ribbonlabels <- paste(lev)
             vrange <- range(intlev)
             if(!is.null(colmap) && !valuesAreColours) 
               col <- colmap(fac)
           },
           character  = {
             x <- eval.im(factor(x))
             lev <- levels(x)
             nvalues <- length(lev)
             trivial <- (nvalues < 2)
             # ensure all factor levels plotted separately
             fac <- factor(lev, levels=lev)
             intlev <- as.integer(fac)
             imagebreaks <- c(intlev - 0.5, max(intlev) + 0.5)
             ribbonvalues <- intlev
             ribbonrange <- range(imagebreaks)
#             ribbonbreaks <- imagebreaks
             ribbonticks <- user.ticks %orifnull% ribbonvalues
             ribbonlabels <- paste(lev)
             vrange <- range(intlev)
             if(!is.null(colmap)) 
               col <- colmap(fac)
           },
           stop(paste("Do not know how to plot image of type", sQuote(xtype)))
           )
  
    ## Compute colour values to be passed to image.default
    if(!is.null(colmap)) {
      ## Explicit colour map object
      colourinfo <- list(breaks=imagebreaks, col=col)
    } else if(!is.null(colfun)) {
      ## Function colfun(n)
      colourinfo <-
        if(is.null(imagebreaks)) list(col=colfun(ncolours)) else
        list(breaks=imagebreaks, col=colfun(length(imagebreaks) - 1L))
    } else if(col.given) {
      ## Colour values
      if(inherits(try(col2rgb(col), silent=TRUE), "try-error"))
        stop("Unable to interpret argument col as colour values")
      if(is.null(imagebreaks)) {
        colourinfo <- list(col=col)
      } else {
        nintervals <- length(imagebreaks) - 1
        colourinfo <- list(breaks=imagebreaks, col=col)
        if(length(col) != nintervals)
          stop(paste("Length of argument", dQuote("col"),
                     paren(paste(length(col))),
                     "does not match the number of distinct values",
                     paren(paste(nintervals))))
      }
    } else stop("Internal error: unable to determine colour values")

    if(spatstat.options("monochrome")) {
      ## transform to grey scale
      colourinfo$col <- to.grey(colourinfo$col)
    }
    
    # colour map to be returned (invisibly)
    i.col <- colourinfo$col
    i.bks <- colourinfo$breaks
    output.colmap <-
      if(is.null(i.col)) NULL else
      if(inherits(i.col, "colourmap")) i.col else
      if(valuesAreColours) colourmap(col=i.col, inputs=i.col) else
      switch(xtype,
             integer=,
             real= {
               if(!is.null(i.bks)) {
                 colourmap(col=i.col, breaks=i.bks)
               } else colourmap(col=i.col, range=vrange, gamma=gamma)
             },
             logical={
               colourmap(col=i.col, inputs=c(FALSE,TRUE))
             },
             character=,
             factor={
               colourmap(col=i.col, inputs=lev)
             },
             NULL)

    ## gamma correction
    soc <- summary(output.colmap)
    if(!is.null(gamma <- soc$gamma) && gamma != 1)
      colourinfo$breaks <- soc$breaks

    ##  ........ decide whether to use rasterImage .........
    
    ## get device capabilities
    ##      (this will start a graphics device if none is active)
    rasterable <- dev.capabilities()$rasterImage
    if(is.null(rasterable)) rasterable <- "no"
    ##
    can.use.raster <-
      switch(rasterable,
             yes=TRUE,
             no=FALSE,
             "non-missing"=!anyNA(x$v),
             FALSE)
    if(is.null(useRaster)) {
      useRaster <- can.use.raster
    } else if(useRaster && !can.use.raster) {
        whinge <- "useRaster=TRUE is not supported by the graphics device"
        if(rasterable == "non-missing")
          whinge <- paste(whinge, "for images with NA values")
        warning(whinge, call.=FALSE)
    } 
    
    ## ........ start plotting .................

    if(!identical(ribbon, TRUE) || trivial) {
      ## no ribbon wanted

      attr(output.colmap, "bbox") <- as.rectangle(x)
      if(!do.plot)
        return(output.colmap)

      ## plot image without ribbon
      image.doit(imagedata=list(x=cellbreaks(x$xcol, x$xstep),
                                y=cellbreaks(x$yrow, x$ystep),
                                z=t(x$v)),
                 W=xbox,
                 workaround=workaround,
                 dotargs,
                 list(useRaster=useRaster, add=add, show.all=show.all),
                 colourinfo,
                 list(xlab = "", ylab = ""),
                 list(asp = 1, main = main, axes=FALSE))
##      if(add && show.all)
##        fakemaintitle(x, main, dotargs)

      do.box.etc(Frame(x), add, dotargs)
      
      return(invisible(output.colmap))
    }
    
    # determine plot region
    bb <- owin(x$xrange, x$yrange)
    Width <- diff(bb$xrange)
    Height <- diff(bb$yrange)
    Size <- max(Width, Height)
    switch(ribside,
           right={
             # ribbon to right of image
             bb.rib <- owin(bb$xrange[2] + c(ribsep, ribsep+ribwid) * Size,
                            bb$yrange)
             rib.iside <- 4
           },
           left={
             # ribbon to left of image
             bb.rib <- owin(bb$xrange[1] - c(ribsep+ribwid, ribsep) * Size,
                            bb$yrange)
             rib.iside <- 2
           },
           top={
             # ribbon above image
             bb.rib <- owin(bb$xrange,
                            bb$yrange[2] + c(ribsep, ribsep+ribwid) * Size)
             rib.iside <- 3
           },
           bottom={
             # ribbon below image
             bb.rib <- owin(bb$xrange,
                            bb$yrange[1] - c(ribsep+ribwid, ribsep) * Size)
             rib.iside <- 1
           })
    bb.all <- boundingbox(bb.rib, bb)

    attr(output.colmap, "bbox") <- bb.all
    attr(output.colmap, "bbox.legend") <- bb.rib
    attr(output.colmap, "side.legend") <- rib.iside
    if(!do.plot)
      return(output.colmap)

    pt <- prepareTitle(main)
    
    if(!add) {
      ## establish coordinate system
      do.call.plotfun(plot.owin,
                      resolve.defaults(list(x=bb.all,
                                            type="n",
                                            main=pt$blank),
                                       dotargs),
                      extrargs=graphicsPars("owin"))
    }
    if(show.all) {
      ## plot title centred over main image area 'bb'
      do.call.plotfun(plot.owin,
                      resolve.defaults(list(x=bb,
                                            type="n",
                                            main=main,
                                            add=TRUE,
                                            show.all=TRUE),
                                       dotargs),
                      extrargs=graphicsPars("owin"))
      main <- ""
    }
    # plot image
    image.doit(imagedata=list(x=cellbreaks(x$xcol, x$xstep),
                              y=cellbreaks(x$yrow, x$ystep),
                              z=t(x$v)),
               W=xbox,
               workaround=workaround,
               list(add=TRUE, show.all=show.all),
               dotargs,
               list(useRaster=useRaster),
               colourinfo,
               list(xlab = "", ylab = ""),
               list(asp = 1, main = main))

##    if(add && show.all)
##      fakemaintitle(bb.all, main, ...)
    
    # box or axes for image
    do.box.etc(bb, add, dotargs)

    # plot ribbon image containing the range of image values
    rib.npixel <- length(ribbonvalues) + 1
    switch(ribside,
           left=,
           right={
             # vertical ribbon
             rib.xcoords <- bb.rib$xrange
             rib.ycoords <- seq(from=bb.rib$yrange[1],
                                to=bb.rib$yrange[2],
                                length.out=rib.npixel)
             rib.z <- matrix(ribbonvalues, ncol=1)
             rib.useRaster <- useRaster
           },
           top=,
           bottom={
             # horizontal ribbon
             rib.ycoords <- bb.rib$yrange
             rib.xcoords <- seq(from=bb.rib$xrange[1],
                                to=bb.rib$xrange[2],
                                length.out=rib.npixel)
             rib.z <- matrix(ribbonvalues, nrow=1)
             # bug workaround
             rib.useRaster <- FALSE 
           })
    image.doit(imagedata=list(x=rib.xcoords,
                              y=rib.ycoords,
                              z=t(rib.z)),
               W=bb.rib,
               workaround=workaround,
               list(add=TRUE,
                    show.all=show.all),
               ribargs,
               list(useRaster=rib.useRaster),
               list(main="", sub=""),
               dotargs,
               colourinfo)
    # box around ribbon?
    resol <- resolve.defaults(ribargs, dotargs)
    if(!identical(resol$box, FALSE))
      plot(as.owin(bb.rib), add=TRUE)
    # scale axis for ribbon image
    ribaxis <- !(identical(resol$axes, FALSE) || identical(resol$ann, FALSE))
    if(ribaxis) {
      ribaxis.iside <- rib.iside
      ## check for user-supplied xlim, ylim with reverse order
      ll <- resolve.defaults(ribargs, dotargs, list(xlim=NULL, ylim=NULL))
      xlimflip <- is.numeric(ll$xlim) && (diff(ll$xlim) < 0)
      ylimflip <- is.numeric(ll$ylim) && (diff(ll$ylim) < 0)
      if(xlimflip) ribaxis.iside <- c(1, 4, 3, 2)[ribaxis.iside] 
      if(ylimflip) ribaxis.iside <- c(3, 2, 1, 4)[ribaxis.iside]
      ##
      axisargs <- list(side=ribaxis.iside, labels=ribbonlabels)
      switch(ribside,
             right={
               scal <- diff(bb.rib$yrange)/diff(ribbonrange)
               at <- bb.rib$yrange[1] + scal * (ribbonticks - ribbonrange[1])
               axisargs <- append(axisargs, list(at=at))
               posargs <- list(pos=bb.rib$xrange[2],
                               yaxp=c(bb.rib$yrange, length(ribbonticks)))
             },
             left={
               scal <- diff(bb.rib$yrange)/diff(ribbonrange)
               at <- bb.rib$yrange[1] + scal * (ribbonticks - ribbonrange[1])
               axisargs <- append(axisargs, list(at=at))
               posargs <- list(pos=bb.rib$xrange[1],
                               yaxp=c(bb.rib$yrange, length(ribbonticks)))
             },
             top={
               scal <- diff(bb.rib$xrange)/diff(ribbonrange)
               at <- bb.rib$xrange[1] + scal * (ribbonticks - ribbonrange[1])
               axisargs <- append(axisargs, list(at=at))
               posargs <- list(pos=bb.rib$yrange[2],
                               xaxp=c(bb.rib$xrange, length(ribbonticks)))
             },
             bottom={
               scal <- diff(bb.rib$xrange)/diff(ribbonrange)
               at <- bb.rib$xrange[1] + scal * (ribbonticks - ribbonrange[1])
               axisargs <- append(axisargs, list(at=at))
               posargs <- list(pos=bb.rib$yrange[1],
                               xaxp=c(bb.rib$xrange, length(ribbonticks)))
             })
      do.call.plotfun(graphics::axis,
                      resolve.defaults(ribargs, axisargs, dotargs, posargs),
                      extrargs=graphicsPars("axis"))
    }
    #
    return(invisible(output.colmap))
  }

  PlotIm
})

invokeColourmapRule <- function(colfun, x, ..., zlim=NULL, colargs=list()) {
  ## utility for handling special functions that generate colour maps
  ## either 
  ##        function(... range) -> colourmap
  ##        function(... inputs) -> colourmap
  stopifnot(is.im(x))
  stopifnot(is.function(colfun))
  colargnames <- names(formals(colfun))
  ## Convert it to a 'colourmap'
  colmap <- NULL
  xtype <- x$type
  if(xtype %in% c("real", "integer") && "range" %in% colargnames) {
    ## function(range) -> colourmap
    vrange <- range(range(x, finite=TRUE), zlim)
    cvals <- try(do.call.matched(colfun,
                                 append(list(range=vrange), colargs)),
                 silent=TRUE)
    if(!inherits(cvals, "try-error")) {
      colmap <- if(inherits(cvals, "colourmap")) cvals else
      if(is.character(cvals)) colourmap(cvals, range=vrange) else NULL
    }
  } else if(xtype != "real" && "inputs" %in% colargnames) {
    ## function(inputs) -> colourmap
    vpossible <- switch(xtype,
                        logical = c(FALSE, TRUE),
                        factor = levels(x),
                        unique(as.matrix(x)))
    if(!is.null(vpossible) && length(vpossible) < 256) {
      cvals <- try(do.call.matched(colfun,
                                   append(list(inputs=vpossible),
                                          colargs)),
                   silent=TRUE)
      if(!inherits(cvals, "try-error")) {
        colmap <- if(inherits(cvals, "colourmap")) cvals else
        if(is.character(cvals))
          colourmap(cvals, inputs=vpossible) else NULL
      }
    }
  }
  return(colmap)
}

########################################################################

image.im <- plot.im

######################################################################

contour.im <- function (x, ..., main, axes=FALSE, add=FALSE,
                        col=par("fg"), 
                        clipwin=NULL, show.all=!add, do.plot=TRUE)
{
  defaultmain <- deparse(substitute(x))
  dotargs <- list(...)
  bb <- Frame(x)
  ## return value
  result <- bb
  attr(result, "bbox") <- bb
  if(!do.plot) return(result)
  ## main title
  sop <- spatstat.options("par.contour")
  if(missing(main)) 
    main <- resolve.1.default(list(main=defaultmain), sop)
  pt <- prepareTitle(main)
  ## plotting parameters
  if(missing(add)) {
    force(add) ## use default in formal arguments, unless overridden
    add <- resolve.1.default(list(add=add), sop)
  }
  if(missing(axes)) {
    force(axes)
    axes <- resolve.1.default(list(axes=axes), sop)
  }
  axes <- axes && !add
  col0 <- if(inherits(col, "colourmap")) par("fg") else col
  ## clip to subset
  if(!is.null(clipwin))
    x <- x[clipwin, drop=FALSE]
  #' start plotting
  if(!add) {
    ## new plot - establish coordinate system
    if(axes && show.all) {
      #' standard plot initialisation in base graphics
      do.call.plotfun(plot.default,
                      resolve.defaults(
                                       list(x = range(x$xcol),
                                            y = range(x$yrow),
                                            type = "n"),
                                       list(...),
                                       list(asp = 1,
                                            xlab = "x",
                                            ylab = "y",
                                            col = col0,
                                            main = main)))
    } else {
      #' plot invisible bounding box
      do.call.plotfun(plot.owin,
                      resolve.defaults(list(x=bb,
                                            type="n",
                                            main=pt$blank),
                                       dotargs),
                      extrargs=graphicsPars("owin"))
    }
  } 
  if(show.all && !axes) {
    ## plot title centred over contour region
    do.call.plotfun(plot.owin,
                    resolve.defaults(list(x=bb,
                                          main=main,
                                          add=TRUE,
                                          show.all=TRUE),
                                     dotargs,
                                     list(col.main=col0)),
                    extrargs=graphicsPars("owin"))
  }
  #' plot contour lines
  if(!inherits(col, "colourmap")) {
    do.call.plotfun(contour.default,
                    resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                     list(add=TRUE, col=col),
                                     list(...)))
  } else {
    clin <- do.call.matched(contourLines,
                            append(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                   list(...)))
    linpar <- graphicsPars("lines")
    for(i in seq_along(clin)) {
      lini <- clin[[i]]
      levi <- lini$level
      coli <- col(levi)
      argi <- resolve.defaults(lini[c("x", "y")],
                               list(...),
                               list(col=coli))
      do.call.matched(lines.default, argi, extrargs=linpar)
    }
  }
  return(invisible(result))
}

