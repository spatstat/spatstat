#
#	plot.owin.S
#
#	The 'plot' method for observation windows (class "owin")
#
#	$Revision: 1.41 $	$Date: 2013/10/09 00:50:18 $
#
#
#

plot.owin <- function(x, main, add=FALSE, ..., box, edge=0.04,
                      hatch=FALSE, angle=45, spacing=diameter(x)/50,
                      invert=FALSE)
{
#
# Function plot.owin.  A method for plot.
#
  if(missing(main))
    main <- short.deparse(substitute(x))
  W <- x
  verifyclass(W, "owin")
  if(missing(box))
    box <- (W$type == "mask")
  else
    stopifnot(is.logical(box) && length(box) == 1)

####
  if(is.expression(main)) 
    nlines <- 1
  else {
    # convert to string and count number of lines
    main <- paste(main)
    if(length(main) > 1)
      main <- paste(main, collapse="\n")
    if(nchar(main) == 0)
      nlines <- 0
    else
      nlines <- length(strsplit(main, "\n")[[1]])
  }
#########        
  xlim <- xr <- W$xrange
  ylim <- yr <- W$yrange

####################################################  
  if(!add) {
    # new plot
    # allow space for main title
    guesslinespace <- 0.08 * diff(yr)
    ylim[2] <- ylim[2] + nlines * guesslinespace
    # set up plot with equal scales
    do.call.matched("plot.default",
                    resolve.defaults(list(x=numeric(0), y=numeric(0),
                                          type="n"),
                                     list(...),
                                     list(xlim=xlim, ylim=ylim,
                                          ann=FALSE, axes=FALSE,
                                          asp=1.0)))
    
    # add title in a reasonable place!
    if(nlines > 0) {
      parval <- resolve.defaults(list(...), par())
      mainheight <- strheight(main, units="user", cex=parval$cex.main)
      gapheight <- (strheight("b\nb", units="user", cex=parval$cex.main)
                    - 2 * strheight("b", units="user", cex=parval$cex.main))
      text(x=mean(xr), y=yr[2] + mainheight + 0.5 * gapheight, labels=main,
           cex=parval$cex.main,
           col=parval$col.main,
           font=parval$font.main)
    }

  }
  
# Draw surrounding box
  if(box)
    do.call.matched("segments",
                    resolve.defaults(
                                     list(x0=xr[c(1,2,2,1)],
                                          y0=yr[c(1,1,2,2)],
                                          x1=xr[c(2,2,1,1)],
                                          y1=yr[c(1,2,2,1)]),
                                     list(...)))
  
# If type = "n", do not plot the window.
    type <- resolve.defaults(list(...), list(type=NULL))$type
    if(!is.null(type) && type == "n")
      return(invisible(NULL))

  
# Draw window

  switch(W$type,
         rectangle = {
           Wpoly <- as.polygonal(W)
           po <- Wpoly$bdry[[1]]
           do.call.matched("polygon",
                           resolve.defaults(list(x=po),
                                            list(...)),
                           extrargs="lwd")
           if(hatch) {
             L <- rlinegrid(angle, spacing, W)
             L <- L[W]
             do.call.matched("plot.psp",
                             resolve.defaults(list(x=L, add=TRUE),
                                              list(...)),
                             extrargs=c("lwd","lty","col"))
           }
         },
         polygonal = {
           p <- W$bdry
           # Determine whether user wants to fill the interior
           col.poly <- resolve.defaults(list(...), list(col=NA))$col
           den.poly <- resolve.defaults(list(...), list(density=NULL))$density
           no.fill  <- is.null(den.poly) &&
                       (is.null(col.poly) || is.na(col.poly))
           # Determine whether we need to triangulate the interior.
           # If it is required to fill the interior,
           # this can be done directly using polygon() provided
           # there are no holes. Otherwise we must triangulate the interior.
           if(no.fill)
             triangulate <- FALSE
           else {
             # Determine whether there are any holes
             holes <- unlist(lapply(p, is.hole.xypolygon))
             triangulate <- any(holes)
           }

           if(!triangulate) {
             # No triangulation required;
             # simply plot the polygons
             for(i in seq_along(p))
               do.call.matched("polygon",
                               resolve.defaults(
                                                list(x=p[[i]]),
                                                list(...)),
                               extrargs="lwd")
           } else {
              # Try using polypath():
             lucy <- names(dev.cur())
             if(!(lucy %in% c("xfig","pictex","X11"))) {
               xx <- unlist(lapply(p, function(a) {c(NA, a$x)}))[-1]
               yy <- unlist(lapply(p, function(a) {c(NA, a$y)}))[-1]
               do.call.matched("polypath",
                               resolve.defaults(list(x=xx,y=yy),
                                                list(border=col.poly),
                                                list(...)))
             } else {
               # decompose window into simply-connected pieces
               broken <- try(break.holes(W))
               if(inherits(broken, "try-error")) {
                 warning("Unable to plot filled polygons")
               } else {
                 # Fill pieces with colour (and draw border in same colour)
                 pp <- broken$bdry
                 for(i in seq_len(length(pp)))
                   do.call.matched("polygon",
                                   resolve.defaults(list(x=pp[[i]],
                                                         border=col.poly),
                                                    list(...)))
               }
             }
             # Now draw polygon boundaries
             for(i in seq_along(p))
               do.call.matched("polygon",
                               resolve.defaults(
                                                list(x=p[[i]]),
                                                list(density=0, col=NA),
                                                list(...)),
                               extrargs="lwd")
           }
           if(hatch) {
             L <- rlinegrid(angle, spacing, W)
             L <- L[W]
             do.call.matched("plot.psp",
                             resolve.defaults(list(x=L, add=TRUE),
                                              list(...)),
                             extrargs=c("lwd","lty","col"))
           }
         },
         mask = {
           # capture 'col' argument and ensure it's at least 2 values
           coldefault <- c(par("bg"), par("fg"))
           col <- resolve.defaults(
                                   list(...),
                                   spatstat.options("par.binary"),
                                   list(col=coldefault)
                                   )$col
           if(length(col) == 1) {
             col <- unique(c(par("bg"), col))
             if(length(col) == 1) 
               col <- c(par("fg"), col)
           }
           # invert colours?
           if(invert)
             col <- rev(col)
           
           do.call.matched("image.default",
                           resolve.defaults(
                           list(x=W$xcol, y=W$yrow, z=t(W$m), add=TRUE),
                           list(col=col),       
                           list(...),
                           spatstat.options("par.binary"),
                           list(zlim=c(FALSE, TRUE))))
           if(hatch)
             warning("Hatching is not implemented for mask windows")
         },
         stop(paste("Don't know how to plot window of type", sQuote(W$type)))
         )

  invisible()
}

break.holes <- function(x, splitby=NULL, depth=0, maxdepth=100) {
  if(is.null(splitby)) {
    # first call: validate x
    stopifnot(is.owin(x))
    splitby <- "x"
  }
  if(depth > maxdepth)
    stop("Unable to divide window into simply-connected pieces")
  p <- x$bdry
  holes <- unlist(lapply(p, is.hole.xypolygon))
  if(!any(holes)) return(x)
  nholes <- sum(holes)
  i <- min(which(holes))
  p.i <- p[[i]]
  b <- as.rectangle(x)
  xr <- b$xrange
  yr <- b$yrange
  switch(splitby,
         x = {
           xsplit <- mean(range(p.i$x))
           left <- c(xr[1], xsplit)
           right <- c(xsplit, xr[2])
           pleft <- intersect.owin(x, owin(left, yr))$bdry
           pright <- intersect.owin(x, owin(right, yr))$bdry
           xnew <- owin(poly=c(pleft, pright), check=FALSE)
           nextsplit <- "y"
         },
         y = {
           ysplit <- mean(range(p.i$y))
           lower <- c(yr[1], ysplit)
           upper <- c(ysplit, yr[2])
           plower <- intersect.owin(x, owin(xr, lower))$bdry
           pupper <- intersect.owin(x, owin(xr, upper))$bdry
           xnew <- owin(poly=c(plower, pupper), check=FALSE)
           nextsplit <- "x"
         })
  # recurse
  xnew <- break.holes(xnew, splitby=nextsplit,
                      depth=depth+1, maxdepth=max(maxdepth, 4*nholes))
  return(xnew)
}
