#'
#'  plot.psp.R
#'
#'  plot method for segment patterns
#'
#'  $Revision: 1.4 $ $Date: 2020/11/17 03:47:24 $

plot.psp <- function(x, ..., main, add=FALSE,
                     show.all=!add, 
                     show.window=show.all,
                     do.plot=TRUE,
                     which.marks=1,
                     style=c("colour", "width", "none"),
                     col=NULL,
                     ribbon=show.all, ribsep=0.15, ribwid=0.05, ribn=1024,
                     scale=NULL, adjust=1,
                     legend=TRUE,
                     leg.side=c("right", "left", "bottom", "top"),
                     leg.sep=0.1,
                     leg.wid=0.1,
                     leg.args=list(),
                     leg.scale=1,
                     negative.args=list(col=2)) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  verifyclass(x, "psp")
  #'
  n <- nsegments(x)
  marx <- marks(x)
  #'
  style <- match.arg(style)
  use.marks <- !is.null(marx) && (n != 0) && (style != "none")
  #'
  if(use.marks && style == "width") {
    #' plot marks as line width
    if(length(dim(marx))) {
      check.1.integer(which.marks)
      marx <- marx[,which.marks]
    }
    values <- as.numeric(marx)
    out <- thickSegments(x, widths=values, ...,
                         add=add,
                         main=main, 
                         do.plot=do.plot,
                         show.all=show.all,
                         show.window=show.window,
                         col=col,
                         negative.args=negative.args,
                         legend=legend,
                         scale=scale, adjust=adjust,
                         leg.side=leg.side,
                         leg.sep=leg.sep,
                         leg.wid=leg.wid,
                         leg.args=leg.args,
                         leg.scale=leg.scale)
    return(invisible(out))
  }
  #' plot marks as colours, if present
  do.ribbon <- identical(ribbon, TRUE) && use.marks
  ##
  ## ....   initialise plot; draw observation window  ......
  owinpars <- setdiff(graphicsPars("owin"), "col")
  if(!do.ribbon) {
    ## window of x only
    bb.all <- as.rectangle(as.owin(x))
    if(do.plot && (!add || show.window)) {
      xwindow <- x$window
      do.call.plotfun(plot.owin, 
                      resolve.defaults(list(x=quote(xwindow),
		                            main=if(show.all) main else "",
                                            add=add,
                                            type = if(show.window) "w" else "n",
                                            show.all=show.all),
                                       list(...)),
                      extrargs=owinpars)
    }
  } else {
    ## enlarged window with room for colour ribbon
    ## x at left, ribbon at right
    bb <- as.rectangle(as.owin(x))
    xwidth <- diff(bb$xrange)
    xheight <- diff(bb$yrange)
    xsize <- max(xwidth, xheight)
    bb.rib <- owin(bb$xrange[2] + c(ribsep, ribsep+ribwid) * xsize,
                   bb$yrange)
    bb.all <- boundingbox(bb.rib, bb)
    if(do.plot) {
      pt <- prepareTitle(main)
      ## establish coordinate system
      if(!add)
      do.call.plotfun(plot.owin,
                      resolve.defaults(list(x=quote(bb.all),
                                            type="n",
                                            main=pt$blank),
                                       list(...)),
                      extrargs=owinpars)
      ## now plot window of x
      ## with title centred on this window
      if(show.window) {
	xwindow <- x$window
        do.call.plotfun(plot.owin, 
                        resolve.defaults(list(x=quote(xwindow),
                                              add=TRUE,
                                              main=main,
                                              show.all=TRUE),
                                         list(...)),
                        extrargs=owinpars)
        ## title done. 
        main <- ""
      }
    }
  }

  # plot segments
  if(n == 0) {
    result <- symbolmap()
    attr(result, "bbox") <- bb.all
    return(invisible(result))
  }
  
  ## determine colours if any
  colmap <- NULL
  if(use.marks) {
    ## use colours
    marx <- as.data.frame(marx)[, which.marks]
    if(is.character(marx) || length(unique(marx)) == 1)
      marx <- factor(marx)
    if(is.null(col)) {
      ## no colour info: use default colour palette
      nc <- if(is.factor(marx)) {
              length(levels(marx))
            } else {
              min(256, length(unique(marx)))
            }
      colfun <- spatstat.options("image.colfun")
      col <- colfun(nc)
    }
    ## determine colour map
    if(inherits(col, "colourmap")) {
      colmap <- colourmap
    } else if(is.colour(col)) {
      ## colour values given; create colour map
      if(is.factor(marx)) {
        lev <- levels(marx)
        colmap <- colourmap(col=col, inputs=factor(lev))
      } else {
        if(!all(is.finite(marx)))
          warning("Some mark values are infinite or NaN or NA")
        colmap <- colourmap(col=col, range=range(marx, finite=TRUE))
      }
    } else stop("Format of argument 'col' is not recognised")
    #' map the mark values to colours
    col <- colmap(marx)
  }
  ## convert to greyscale?
  if(spatstat.options("monochrome")) {
    col <- to.grey(col)
    colmap <- to.grey(colmap)
  }
  if(do.plot) {
    ## plot segments
    do.call.plotfun(segments,
                    resolve.defaults(as.list(x$ends),
                                     list(...),
                                     list(col=col),
                                     .MatchNull=FALSE, .StripNull=TRUE),
                    extrargs=names(par()))
    ## plot ribbon
    if(do.ribbon) 
      plot(colmap, vertical=TRUE, add=TRUE,
           xlim=bb.rib$xrange, ylim=bb.rib$yrange)
  }
  
  # return colour map
  result <- colmap %orifnull% colourmap()
  attr(result, "bbox") <- bb.all
  return(invisible(result))
}

thickSegments <- local({

  ## plot segment pattern with variable widths
  
  thickSegments <- function(x, widths, ...,
                            add=FALSE,
                            main="",
                            do.plot=TRUE,
                            show.all=!add, 
                            show.window=show.all,
                            scale=NULL, adjust=1, 
                            negative.args=list(col=2),
                            legend=TRUE,
                            leg.side=c("right", "left", "bottom", "top"),
                            leg.sep=0.1,
                            leg.wid=0.1,
                            leg.args=list(),
                            leg.scale=1,
                            zlim,
                            box=FALSE) {
    leg.side <- match.arg(leg.side)
    check.1.real(leg.scale)
    check.1.real(adjust)

    if(missing(zlim) || is.null(zlim)) {
      zlim <- NULL
      zliminfo <- list()
    } else {
      check.range(zlim)
      stopifnot(all(is.finite(zlim)))
      zliminfo <- list(zlim=zlim)
    }
    
    W <- Window(x)
    stopifnot(is.numeric(widths))
    #' convert non-finite widths to zero width
    widths[!is.finite(widths)] <- 0
    #' plan layout
    if(legend) {
      #' use layout procedure in plot.im
      px <- pixellate(x)
      z <- do.call(plot.im,
                   resolve.defaults(list(quote(px), 
					do.plot=FALSE, ribbon=TRUE),
                                    list(...),
                                    list(ribside  = leg.side,
                                         ribsep   = leg.sep,
                                         ribwid   = leg.wid,
                                         ribargs  = leg.args,
                                         ribscale = leg.scale),
                                    zliminfo,
                                    list(main=main, valuesAreColours=FALSE)))
      bb.all <- attr(z, "bbox")
      bb.leg <- attr(z, "bbox.legend")
    } else {
      bb.all <- Frame(W)
      bb.leg <- NULL
    }
    legend <- !is.null(bb.leg)
    if(legend) {
      #' expand plot region to accommodate text annotation in legend
      if(leg.side %in% c("left", "right")) {
        delta <- 2 * sidelengths(bb.leg)[1L]
        xmargin <- if(leg.side == "right") c(0, delta) else c(delta, 0)
        bb.all <- grow.rectangle(bb.all, xmargin=xmargin)
      }
    }
    #' initialise plot
    bb <- do.call.matched(plot.owin,
                          resolve.defaults(list(x=quote(bb.all), type="n"),
                                           list(...), list(main=main)),
                          extrargs="type")
    if(box)
      plot(Frame(W), add=TRUE)
    #' resolve graphics parameters for polygons
    names(negative.args) <- paste0(names(negative.args), ".neg")
    grafpar <- resolve.defaults(negative.args,
                                list(...),
                               list(col=1),
                                .MatchNull=FALSE)
    #' rescale width values to a plottable range
    if(is.null(zlim)) zlim <- range(widths, finite=TRUE)
    vr <- range(0, zlim)
    if(is.null(scale)) {
      maxsize <- mean(distmap(x))/2
      scale <- maxsize/max(abs(vr))
    } else check.1.real(scale)
    phys.scale <- adjust * scale
    halfwidths <- phys.scale * widths/2
    #' plot each segment
    thetaperp <- angles.psp(x) + pi/2
    ends <- as.matrix(unclass(x)$ends)
    for(i in seq_len(nobjects(x))) {
      xx <- ends[i, c(1L,3L)]
      yy <- ends[i, c(2L,4L)]
      drawseg(xx, yy, rep.int(halfwidths[i], 2L), thetaperp[i], grafpar)
    }
    result <- phys.scale 
    attr(result, "bbox") <- bb
    if(legend) {
      attr(result, "bbox.legend") <- bb.leg
      plotWidthMap(bb.leg     = bb.leg,
                   zlim       = zlim,
                   phys.scale = phys.scale,
                   leg.scale  = leg.scale,
                   leg.side   = leg.side,
                   leg.args   = leg.args,
                   grafpar    = grafpar)
    }
    return(invisible(result))
  }

  drawseg <- function(xx, yy, vv, ang, pars) {
    ## draw polygon around segment
    sgn <- sign(mean(vv))
    xx <- c(xx, rev(xx))
    yy <- c(yy, rev(yy))
    vv <- c(vv, -rev(vv))
    xx <- xx + cos(ang) * vv
    yy <- yy + sin(ang) * vv
    drawSignedPoly(xx, yy, pars, sgn)
    invisible(NULL)
  }

  thickSegments

})

drawSignedPoly <- local({
  
  ## internal function to plot line segments for style="width"
  ## with sign-dependent colours, etc

  pNames <- c("density", "angle", "border", "col", "lty")
  posnames <- paste(pNames, ".pos", sep="")
  negnames <- paste(pNames, ".neg", sep="")
  
  redub <- function(from, to, x) {
    #' rename entry x$from to x$to
    m <- match(from, names(x))
    if(any(ok <- !is.na(m))) 
      names(x)[m[ok]] <- to[ok]
    return(resolve.defaults(x))
  }
  
  drawSignedPoly <- function(x, y, pars, sgn) {
    #' plot polygon using parameters appropriate to "sign"
    if(sgn >= 0) {
      pars <- redub(posnames, pNames, pars)
    } else {
      pars <- redub(negnames, pNames, pars)
    }
    pars <- pars[names(pars) %in% pNames]
    if(is.null(pars$border)) pars$border <- pars$col
    do.call(polygon, append(list(x=x, y=y), pars))
    invisible(NULL)
  }

  drawSignedPoly
})

## internal function to plot the map of pixel values to line widths

plotWidthMap <- function(bb.leg, zlim, phys.scale,
                         leg.scale, leg.side,
                         leg.args, grafpar) {
  ## get graphical arguments
  grafpar <- resolve.defaults(leg.args, grafpar)
  ## set up scale of typical pixel values
  gvals <- leg.args$at %orifnull% prettyinside(zlim)
  ## corresponding widths
  wvals <- phys.scale * gvals
  ## glyph positions
  ng <- length(gvals)
  xr <- bb.leg$xrange
  yr <- bb.leg$yrange
  switch(leg.side,
         right = ,
         left = {
           y <- seq(yr[1], yr[2], length.out=ng+1L)
           y <- (y[-1L] + y[-(ng+1L)])/2
           for(j in 1:ng) {
             xx <- xr[c(1L,2L,2L,1L)]
             yy <- (y[j] + c(-1,1) * wvals[j]/2)[c(1L,1L,2L,2L)]
             drawSignedPoly(x = xx, y = yy, grafpar, sign(wvals[j]))
           }
         },
         bottom = ,
         top = {
           x <- seq(xr[1], xr[2], length.out=ng+1L)
           x <- (x[-1L] + x[-(ng+1L)])/2
           for(j in 1:ng) {
             xx <- (x[j] + c(-1,1) * wvals[j]/2)[c(1L,1L,2L,2L)]
             yy <- yr[c(1L,2L,2L,1L)]
             drawSignedPoly(x = xx, y = yy, grafpar, sign(wvals[j]))
           }
         })
  ## add text labels
  glabs <- signif(leg.scale * gvals, 2)
  textpos <- switch(leg.side,
                    right  = list(x=xr[2], y=y,     pos=4),
                    left   = list(x=xr[1], y=y,     pos=2),
                    bottom = list(x=x,     y=yr[1], pos=1),
                    top    = list(x=x,     y=yr[2], pos=3))
  textargs <- resolve.defaults(textpos,
                               leg.args,
                               list(labels=glabs))
  do.call.matched(text, textargs, extrargs=graphicsPars("text"))
  return(invisible(NULL))
}
