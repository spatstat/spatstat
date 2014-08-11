#
#   plot.im.R
#
#  $Revision: 1.67 $   $Date: 2013/07/17 05:51:40 $
#
#  Plotting code for pixel images
#
#  plot.im
#  image.im
#  contour.im
#  persp.im
#
###########################################################################

plot.im <- local({

  # recognised additional arguments to image.default() and axis()

  plotparams <- imageparams <- c("main", "asp", "sub", "axes", "ann",
                   "cex", "font", 
                   "cex.axis", "cex.lab", "cex.main", "cex.sub",
                   "col.axis", "col.lab", "col.main", "col.sub",
                   "font.axis", "font.lab", "font.main", "font.sub")

  axisparams <- c("cex", 
                  "cex.axis", "cex.lab",
                  "col.axis", "col.lab",
                  "font.axis", "font.lab")
  
  # auxiliary functions

  image.doit <- function(..., extrargs=imageparams) {
    do.call.matched("image.default",
                    resolve.defaults(...),
                    extrargs=extrargs)
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
    ok <- (x > 0)
    y[ok] <- log10(x[ok])
    return(y)
  }
  
  # main function
  PlotIm <- function(x, ...,
                     col=NULL, valuesAreColours=NULL, log=FALSE,
                     ribbon=TRUE, ribside=c("right", "left", "bottom", "top"),
                     ribsep=0.15, ribwid=0.05, ribn=1024,
                     ribscale=1, ribargs=list(), colargs=list()) {
    main <- short.deparse(substitute(x))
    verifyclass(x, "im")
    dotargs <- list(...)
    stopifnot(is.list(ribargs))
    ribside <- match.arg(ribside)
    zlim <- dotargs$zlim

    x <- repair.image.xycoords(x)

    xtype <- x$type

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
        } else if(!is.null(col)) {
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
    } else if(!is.null(col)) {
      # argument 'col' controls colours
      valuesAreColours <- FALSE
    } else {
      # default : determine whether pixel values are colours
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
      rx <- range(x)
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
    ribbonvalues <- ribbonbreaks <- NULL

    # predetermined colour map?
    # (i.e. mapping from values to colours)
    colmap <- if(inherits(col, "colourmap")) col else NULL

    # colour map determined by a function?
    if(is.null(colmap) && (is.null(col) || is.function(col))) {
      colfun <- if(is.function(col)) col else spatstat.options("image.colfun")
      colargnames <- names(formals(colfun))
      if("range" %in% colargnames && xtype %in% c("real", "integer")) {
        # continuous 
        vrange <- range(range(x), zlim)
        cvals <- try(do.call.matched(colfun,
                                     append(list(range=vrange), colargs)),
                     silent=TRUE)
        if(!inherits(cvals, "try-error")) {
          colmap <- if(inherits(cvals, "colourmap")) cvals else
            if(is.character(cvals)) colourmap(cvals, range=vrange) else NULL
        }
      } else if("inputs" %in% colargnames && xtype != "real") {
        # discrete
        vpossible <- switch(xtype,
                            logical = c(FALSE, TRUE),
                            factor = levels(x),
                            unique(as.matrix(x)))
        if(!is.null(vpossible) && length(vpossible) < 256) {
          cvals <- try(do.call.matched(colfun,
                                       append(list(inputs=vpossible), colargs)),
                       silent=TRUE)
          if(!inherits(cvals, "try-error")) {
            colmap <- if(inherits(cvals, "colourmap")) cvals else
            if(is.character(cvals)) colourmap(cvals, inputs=vpossible) else NULL
          }
        }
      }
    }
       
    switch(xtype,
           real    = {
             vrange <- range(x)
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
               nominalmarks <- axisTicks(nominalrange, log=do.log)
               ribbonticks <- Log(nominalmarks/ribscale)
               ribbonlabels <- paste(nominalmarks)
             }
           },
           integer = {
             values <- as.vector(x$v)
             values <- values[!is.na(values)]
             uv <- unique(values)
             vrange <- range(uv)
             vrange <- range(zlim, vrange)
             nvalues <- length(uv)
             trivial <- (nvalues < 2)
             if(!trivial){
               nominalrange <- Log(ribscale * Exp(vrange))
               nominalmarks <- axisTicks(nominalrange, log=do.log)
               nominalmarks <- nominalmarks[nominalmarks %% 1 == 0]
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
             ribbonbreaks <- imagebreaks
             ribbonticks <- ribbonvalues
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
             ribbonbreaks <- imagebreaks
             ribbonticks <- ribbonvalues
             ribbonlabels <- paste(lev)
             vrange <- range(intlev)
             if(!is.null(colmap)) 
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
             ribbonbreaks <- imagebreaks
             ribbonticks <- ribbonvalues
             ribbonlabels <- paste(lev)
             vrange <- range(intlev)
             if(!is.null(colmap)) 
               col <- colmap(fac)
           },
           stop(paste("Do not know how to plot image of type", sQuote(xtype)))
           )
  
    # determine colour map
    if(!is.null(colmap)) {
      # explicit colour map object
      colourinfo <- list(breaks=imagebreaks, col=col)
    } else {
      # compile colour information
      # start with default colour values
      colfun <- spatstat.options("image.colfun")
      # we already know that colfun(n) works
      colourinfo <- if(!is.null(imagebreaks)) {
        list(breaks=imagebreaks, col=colfun(length(imagebreaks)-1))
      } else list(col=colfun(256))
      if(!is.null(col) && is.character(col)) {
        # overwrite colour info with user-specified colour values
        colourinfo$col <- col
        if(!is.null(colourinfo$breaks)) {
        # check consistency
          nvalues <- length(colourinfo$breaks) - 1
          if(length(col) != nvalues)
            stop(paste("Length of argument", dQuote("col"),
                       paren(paste(length(col))),
                       "does not match the number of distinct values",
                       paren(paste(nvalues))))
        }
      }
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
               } else colourmap(col=i.col, range=vrange)
             },
             logical={
               colourmap(col=i.col, inputs=c(FALSE,TRUE))
             },
             character=,
             factor={
               colourmap(col=i.col, inputs=lev)
             },
             NULL)
    
    # ........ secret exit used by plot.listof
    
    preponly <- resolve.defaults(dotargs, list(preponly=FALSE))$preponly
    if(preponly) return(output.colmap)

    # ........ start plotting .................

    add <- resolve.defaults(dotargs, list(add=FALSE))$add

    if(!identical(ribbon, TRUE)
       || identical(add, TRUE)
       || trivial)
      {
        # plot image without ribbon
        image.doit(list(x=cellbreaks(x$xcol, x$xstep),
                        y=cellbreaks(x$yrow, x$ystep),
                        z=t(x$v)),
                   dotargs,
                   list(useRaster=TRUE),
                   colourinfo,
                   list(xlab = "", ylab = ""),
                   list(asp = 1, main = main, axes=FALSE)
                   )
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
    bb.all <- bounding.box(bb.rib, bb)
    # establish coordinate system
    do.call.matched("plot.default",
                    resolve.defaults(list(x=0, y=0,  type="n",
                                          axes=FALSE, asp=1,
                                          xlim=bb.all$xrange,
                                          ylim=bb.all$yrange),
                                     dotargs,
                                     list(main=main, xlab="", ylab="")),
                    extrargs=plotparams)
    # plot image
    image.doit(list(x=cellbreaks(x$xcol, x$xstep),
                    y=cellbreaks(x$yrow, x$ystep),
                    z=t(x$v)),
               list(add=TRUE),
               dotargs,
               list(useRaster=TRUE),
               colourinfo,
               list(xlab = "", ylab = ""),
               list(asp = 1, main = main))
    # axes for image
    imax <- identical(dotargs$axes, TRUE)
    imbox <- !identical(dotargs$box, FALSE)
    if(imbox)
      rect(x$xrange[1], x$yrange[1], x$xrange[2], x$yrange[2])
    if(imax) {
      px <- pretty(bb$xrange)
      py <- pretty(bb$yrange)
      do.call.matched("axis",
                      resolve.defaults(
                                       list(side=1, at=px), 
                                       dotargs,
                                       list(pos=bb$yrange[1])),
                      extrargs=axisparams)
      do.call.matched("axis",
                      resolve.defaults(
                                       list(side=2, at=py), 
                                       dotargs,
                                       list(pos=bb$xrange[1])),
                      extrargs=axisparams)
    }
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
             rib.useRaster <- TRUE
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
    image.doit(list(x=rib.xcoords, y=rib.ycoords,
                    z=t(rib.z),
                    add=TRUE),
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
      axisargs <- list(side=rib.iside, labels=ribbonlabels)
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
      do.call.matched("axis",
                      resolve.defaults(axisargs,
                                       ribargs, dotargs,
                                       posargs),
                      extrargs=axisparams)
    }
    #
    return(invisible(output.colmap))
  }

  PlotIm
})



########################################################################

image.im <- plot.im

########################################################################

persp.im <- function(x, ..., colmap=NULL) {
  xname <- deparse(substitute(x))
  xinfo <- summary(x)
  if(xinfo$type == "factor")
    stop("Perspective plot is inappropriate for factor-valued image")
  pop <- spatstat.options("par.persp")
  # check for common error
  if(!is.null(col <- list(...)$col) && !is.matrix(col))
    warning("Argument col is not a matrix. Did you mean colmap?")
  # colour map?
  if(is.null(colmap)) {
    colinfo <- list(col=NULL)
  } else if(inherits(colmap, "colourmap")) {
    # colour map object
    # apply colour function to image data
    colval <- eval.im(colmap(x))
    colval <- t(as.matrix(colval))
    # strip one row and column for input to persp.default
    colval <- colval[-1, -1]
    # replace NA by arbitrary value
    isna <- is.na(colval)
    if(any(isna)) {
      stuff <- attr(colmap, "stuff")
      colvalues <- stuff$outputs
      colval[isna] <- colvalues[1]
    }
    # pass colour matrix (and suppress lines)
    colinfo <- list(col=colval, border=NA)
  } else {
    # interpret 'colmap' as colour map
    if(is.list(colmap) && all(c("breaks", "col") %in% names(colmap))) {
      breaks <- colmap$breaks
      colvalues <- colmap$col
    } else if(is.vector(colmap)) {
      colvalues <- colmap
      breaks <- quantile(x, seq(from=0,to=1,length.out=length(colvalues)+1))
      if(!all(ok <- !duplicated(breaks))) {
        breaks <- breaks[ok]
        colvalues <- colvalues[ok[-1]]
      }
    } else warning("Unrecognised format for colour map")
    # apply colour map to image values
    colid <- cut.im(x, breaks=breaks, include.lowest=TRUE)
    colval <- eval.im(colvalues[unclass(colid)])
    colval <- t(as.matrix(colval))
    nr <- nrow(colval)
    nc <- ncol(colval)
    colval <- colval[-1, -1]
    colval[is.na(colval)] <- colvalues[1]
    # pass colour matrix (and suppress lines)
    colinfo <- list(col=colval, border=NA)
  }

  # get reasonable z scale while fixing x:y aspect ratio
  if(xinfo$type %in% c("integer", "real")) {
    zrange <- xinfo$range
    if(diff(zrange) > 0) {
      xbox <- as.rectangle(x)
      zscale <- 0.5 * mean(diff(xbox$xrange), diff(xbox$yrange))/diff(zrange)
      zlim <- zrange
    } else {
      zscale <- NULL
      zlim <- c(0,2) * xinfo$mean
    }
  } else 
    zscale <- zlim <- NULL

  do.call.matched("persp",
                  resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                   list(...),
                                   pop,
                                   colinfo,
                                   list(xlab="x", ylab="y", zlab=xname),
                                   list(scale=FALSE, expand=zscale, zlim=zlim),
                                   list(main=xname),
                                   .StripNull=TRUE),
                  funargs=.Spatstat.persp.args)
}

.Spatstat.persp.args <- list("x", "y", "z",
                             "xlim", "ylim", "zlim",
                             "xlab", "ylab", "zlab",
                             "main", "sub",
                             "theta", "phi", "r", "d", "scale",
                             "expand", "col", "border",
                             "ltheta", "lphi", "shade", "box",
                             "axes", "nticks", "ticktype")

######################################################################

contour.im <- function (x, ..., main, axes=TRUE, add=FALSE)
{
  defaultmain <- deparse(substitute(x))
  sop <- spatstat.options("par.contour")
  if(missing(main)) 
    main <- resolve.defaults(sop, list(main=defaultmain))$main
  if(missing(add))
    add <- resolve.defaults(sop, list(add=FALSE))$add
  if(missing(axes))
     axes <- resolve.defaults(sop, list(axes=TRUE))$axes
  if(!add) {
    # new plot
    if(axes) # with axes
      do.call.matched("plot.default",
                      resolve.defaults(
                                       list(x = range(x$xcol),
                                            y = range(x$yrow),
                                            type = "n"),
                                       list(...),
                                       list(asp = 1, xlab = "x",
                                            ylab = "y", main = main)))
    else { # box without axes
      rec <- owin(x$xrange, x$yrange)
      do.call.matched("plot.owin",
                      resolve.defaults(list(x=rec),
                                       list(...),
                                       list(main=main)))
    }
  }
  do.call.matched("contour.default",
                  resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                   list(add=TRUE),
                                   list(...)))
  return(invisible(NULL))
}

