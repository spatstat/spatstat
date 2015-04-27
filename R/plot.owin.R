#
#	plot.owin.S
#
#	The 'plot' method for observation windows (class "owin")
#
#	$Revision: 1.55 $	$Date: 2015/04/27 07:24:25 $
#
#
#

plot.owin <- function(x, main, add=FALSE, ..., box, edge=0.04,
                      type = c("w", "n"), show.all=!add,
                      hatch=FALSE,
                      hatchargs=list(), 
                      invert=FALSE, do.plot=TRUE,
                      claim.title.space=FALSE) 
{
#
# Function plot.owin.  A method for plot.
#
  if(missing(main))
    main <- short.deparse(substitute(x))

  W <- x
  verifyclass(W, "owin")
  if(!do.plot) 
    return(invisible(as.rectangle(W)))
  
  type <- match.arg(type)

  if(missing(box) || is.null(box)) {
    box <- is.mask(W) && show.all
  } else stopifnot(is.logical(box) && length(box) == 1)

####
  pt <- prepareTitle(main)
  main <- pt$main
  nlines <- pt$nlines
#########        
  xlim <- xr <- W$xrange
  ylim <- yr <- W$yrange

####################################################

  ## graphics parameters that can be overridden by user
  gparam <- resolve.defaults(list(...), par())
  ## character expansion factors
  ##     main title size = 'cex.main' * par(cex.main) * par(cex)
  ## user's graphics expansion factor (*multiplies* par)
  cex.main.user <- resolve.1.default(list(cex.main=1), list(...))
  ## size of main title as multiple of par('cex')
  cex.main.rela <- cex.main.user * par('cex.main') 
  ## absolute size
  cex.main.absol <- cex.main.rela * par('cex')
    
  if(!add) {
    ## new plot
    if(claim.title.space && nlines > 0) {
      ## allow space for main title (only in multi-panel plots)
      guesslinespace <- 0.07 * sqrt(diff(xr)^2 + diff(yr)^2) * cex.main.absol
      added <- (nlines + 1) * guesslinespace
      ylim[2] <- ylim[2] + added
    }
    ## set up plot with equal scales
    do.call.plotfun("plot.default",
                    resolve.defaults(list(x=numeric(0), y=numeric(0),
                                          type="n"),
                                     list(...),
                                     list(xlim=xlim, ylim=ylim,
                                          ann=FALSE, axes=FALSE,
                                          asp=1.0,
                                          xaxs="i", yaxs="i")))
  }
  if(show.all && nlines > 0) {
    ## add title 
    if(claim.title.space) {
      mainheight <- sum(strheight(main, units="user", cex=cex.main.rela))
      gapheight <- (strheight("b\nb", units="user", cex=cex.main.rela)
                    - 2 * strheight("b", units="user", cex=cex.main.rela))
      if(nlines > 1 && !is.expression(main))
        main <- paste(main, collapse="\n")
      text(x=mean(xr), y=yr[2] + mainheight + 0.5 * gapheight, labels=main,
           cex=cex.main.rela,
           col=gparam$col.main,
           font=gparam$font.main)
    } else {
      title(main=main,
            cex=cex.main.rela,
            col=gparam$col.main,
            font=gparam$font.main)
    }
  }
  
# Draw surrounding box
  if(box)
    do.call.plotfun("segments",
                    resolve.defaults(
                                     list(x0=xr[c(1,2,2,1)],
                                          y0=yr[c(1,1,2,2)],
                                          x1=xr[c(2,2,1,1)],
                                          y1=yr[c(1,2,2,1)]),
                                     list(...)))
  
# If type = "n", do not plot the window.
    if(type == "n")
      return(invisible(as.rectangle(W)))

  
# Draw window

  switch(W$type,
         rectangle = {
           Wpoly <- as.polygonal(W)
           po <- Wpoly$bdry[[1]]
           do.call.plotfun("polygon",
                           resolve.defaults(list(x=po),
                                            list(...)),
                           extrargs="lwd")
           if(hatch)
             do.call("add.texture", append(list(W=W), hatchargs))
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
               do.call.plotfun("polygon",
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
               do.call.plotfun("polypath",
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
                   do.call.plotfun("polygon",
                                   resolve.defaults(list(x=pp[[i]],
                                                         border=col.poly),
                                                    list(...)))
               }
             }
             # Now draw polygon boundaries
             for(i in seq_along(p))
               do.call.plotfun("polygon",
                               resolve.defaults(
                                                list(x=p[[i]]),
                                                list(density=0, col=NA),
                                                list(...)),
                               extrargs="lwd")
           }
           if(hatch)
             do.call("add.texture", append(list(W=W), hatchargs))
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
           ## invert colours?
           if(invert)
             col <- rev(col)
           ## convert to greyscale?
           if(spatstat.options("monochrome"))
             col <- to.grey(col)
           
           do.call.matched("image.default",
                           resolve.defaults(
                           list(x=W$xcol, y=W$yrow, z=t(W$m), add=TRUE),
                           list(col=col),       
                           list(...),
                           spatstat.options("par.binary"),
                           list(zlim=c(FALSE, TRUE))))
           if(hatch)
             do.call("add.texture", append(list(W=W), hatchargs))
         },
         stop(paste("Don't know how to plot window of type", sQuote(W$type)))
         )
  return(invisible(as.rectangle(W)))
}

break.holes <- local({

  insect <- function(A, Box) {
    ## efficient version of intersect.owin which doesn't 'fix' the polygons
    a <- lapply(A$bdry, reverse.xypolygon)
    b <- lapply(as.polygonal(Box)$bdry, reverse.xypolygon)
    ab <- polyclip::polyclip(a, b, "intersection",
                             fillA="nonzero", fillB="nonzero")
    if(length(ab)==0)
      return(emptywindow(Box))
    # ensure correct polarity
    totarea <- sum(unlist(lapply(ab, Area.xypolygon)))
    if(totarea < 0)
      ab <- lapply(ab, reverse.xypolygon)
    AB <- owin(Box$xrange, Box$yrange,
               poly=ab, check=FALSE, strict=FALSE, fix=FALSE,
               unitname=unitname(A))
    return(AB)
  }

  break.holes <- function(x, splitby=NULL, depth=0, maxdepth=100) {
    if(is.null(splitby)) {
      ## first call: validate x
      stopifnot(is.owin(x))
      splitby <- "x"
    }
    if(depth > maxdepth)
      stop("Unable to divide window into simply-connected pieces")
    p <- x$bdry
    holes <- unlist(lapply(p, is.hole.xypolygon))
    if(!any(holes)) return(x)
    nholes <- sum(holes)
    maxdepth <- max(maxdepth, 4 * nholes)
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
             xleft <- insect(x, owin(left, yr))
             xright <- insect(x, owin(right, yr))
             ## recurse
             xleft <- break.holes(xleft, splitby="y",
                                  depth=depth+1, maxdepth=maxdepth)
             xright <- break.holes(xright, splitby="y",
                                  depth=depth+1, maxdepth=maxdepth)
             ## recombine (without fusing polygons again!)
             result <- owin(xr, yr, poly=c(xleft$bdry, xright$bdry),
                            check=FALSE, strict=FALSE, fix=FALSE)
           },
           y = {
             ysplit <- mean(range(p.i$y))
             lower <- c(yr[1], ysplit)
             upper <- c(ysplit, yr[2])
             xlower <- insect(x, owin(xr, lower))
             xupper <- insect(x, owin(xr, upper))
             ## recurse
             xlower <- break.holes(xlower, splitby="x",
                                   depth=depth+1, maxdepth=maxdepth)
             xupper <- break.holes(xupper, splitby="x",
                                   depth=depth+1, maxdepth=maxdepth)
             ## recombine (without fusing polygons again!)
             result <- owin(xr, yr, poly=c(xlower$bdry, xupper$bdry),
                            check=FALSE, strict=FALSE, fix=FALSE)
           })
    return(result)
  }

  break.holes
})

