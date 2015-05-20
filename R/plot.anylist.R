##
##  plot.anylist.R
##
##  Plotting functions for 'solist', 'anylist', 'imlist'
##       and legacy class 'listof'
##
##  $Revision: 1.15 $ $Date: 2015/05/20 14:41:37 $
##

plot.anylist <- plot.solist <- plot.listof <-
  local({

  ## auxiliary functions
  extraplot <- function(nnn, x, ..., add=FALSE, extrargs=list(),
                        panel.args=NULL, plotcommand="plot") {
    argh <- list(...)
    if(is.ppp(x) && identical(plotcommand,"plot"))
      argh <- c(argh, list(multiplot=FALSE))
    if(!is.null(panel.args)) {
      xtra <- if(is.function(panel.args)) panel.args(nnn) else panel.args
      if(!is.list(xtra))
        stop(paste0("panel.args",
                    if(is.function(panel.args)) "(i)" else "",
                    " should be a list"))
      argh <- resolve.defaults(xtra, argh)
    }
    if(length(extrargs) > 0)
      argh <- resolve.defaults(argh, extrargs)
    ## some plot commands don't recognise 'add'
    if(add)
      argh <- append(argh, list(add=TRUE))
    do.call(plotcommand, append(list(x=x), argh))
  }

  exec.or.plot <- function(cmd, i, xi, ..., extrargs=list(), add=FALSE) {
    if(is.null(cmd)) return(NULL)
    argh <- resolve.defaults(list(...),
                             extrargs,
                             ## some plot commands don't recognise 'add' 
                             if(add) list(add=TRUE) else NULL,
                             if(is.ppp(cmd)) list(multiplot=FALSE) else NULL)
    if(is.function(cmd)) {
      do.call(cmd, resolve.defaults(list(i, xi), argh))
    } else {
      do.call(plot, resolve.defaults(list(cmd), argh))
    }
  }

  exec.or.plotshift <- function(cmd, i, xi, ..., vec=vec,
                                extrargs=list(), add=FALSE) {
    if(is.null(cmd)) return(NULL)
    argh <- resolve.defaults(list(...),
                             extrargs,
                             ## some plot commands don't recognise 'add' 
                             if(add) list(add=TRUE) else NULL,
                             if(is.ppp(cmd)) list(multiplot=FALSE) else NULL)
    if(is.function(cmd)) {
      do.call(cmd, resolve.defaults(list(i, xi), argh))
    } else {
      cmd <- shift(cmd, vec)
      do.call(plot, resolve.defaults(list(cmd), argh))
    }
  }

  ## bounding box, including ribbon for images, legend for point patterns
  getplotbox <- function(x, ..., do.plot, plotcommand="plot") {
    if(inherits(x, c("im", "ppp", "psp", "msr", "layered", "tess"))) {
      if(identical(plotcommand, "plot")) {
        y <- if(is.ppp(x)) plot(x, ..., multiplot=FALSE, do.plot=FALSE) else 
                           plot(x, ..., do.plot=FALSE)      
        return(as.owin(y))
      } else if(identical(plotcommand, "contour")) {
        y <- contour(x, ..., do.plot=FALSE)      
        return(as.owin(y))
      } else {
        plc <- plotcommand
        if(is.character(plc)) plc <- get(plc)
        if(!is.function(plc)) stop("Unrecognised plot function")
        if("do.plot" %in% names(args(plc)))
          return(as.owin(do.call(plc, list(x=x, ..., do.plot=FALSE))))
        return(as.rectangle(x))
      }
    }
    return(try(as.rectangle(x), silent=TRUE))
  }

  # calculate bounding boxes for each panel using intended arguments!
  getPlotBoxes <- function(xlist, ..., panel.args=NULL, extrargs=list()) {
    userargs <- list(...)
    n <- length(xlist)
    result <- vector(length=n, mode="list")
    for(i in seq_len(n)) {
      pai <- if(is.function(panel.args)) panel.args(i) else list()
      argh <- resolve.defaults(pai, userargs, extrargs)
      result[[i]] <- do.call(getplotbox, append(list(x=xlist[[i]]), argh))
    }
    return(result)
  }
    
  is.shiftable <- function(x) {
    if(is.null(x)) return(TRUE)
    if(is.function(x)) return(FALSE)
    y <- try(as.rectangle(x), silent=TRUE)
    return(!inherits(y, "try-error"))
  }

  plot.anylist <- function(x, ..., main, arrange=TRUE,
                            nrows=NULL, ncols=NULL,
                            main.panel=NULL,
                            mar.panel=c(2,1,1,2),
                            hsep = 0,
                            vsep = 0,
                            panel.begin=NULL,
                            panel.end=NULL,
                            panel.args=NULL,
                            plotcommand="plot",
                            adorn.left=NULL,
                            adorn.right=NULL,
                            adorn.top=NULL,
                            adorn.bottom=NULL,
                            adorn.size=0.2,
                            equal.scales=FALSE,
                            halign=FALSE, valign=FALSE) {
    xname <- short.deparse(substitute(x))

    isSo <- inherits(x, "solist")
    isIm <- inherits(x, "imlist") || (isSo && all(unlist(lapply(x, is.im))))
    
    ## `boomerang despatch'
    cl <- match.call()
    if(missing(plotcommand) && isIm) {
      cl[[1]] <- as.name("image.imlist")
      parenv <- sys.parent()
      return(invisible(eval(cl, envir=parenv)))
    }

    if(isSo) {
      allfv <- somefv <- FALSE
    } else {
      isfv <- unlist(lapply(x, is.fv))
      allfv <- all(isfv)
      somefv <- any(isfv)
    }
    
    ## panel margins
    if(!missing(mar.panel)) {
      nm <- length(mar.panel)
      if(nm == 1) mar.panel <- rep(mar.panel, 4) else
      if(nm == 2) mar.panel <- rep(mar.panel, 2) else
      if(nm != 4) stop("mar.panel should have length 1, 2 or 4")
    } else if(somefv) {
      ## change default
      mar.panel <- 0.25+c(4,4,2,2)
    }
    
    n <- length(x)
    names(x) <- good.names(names(x), "Component_", 1:n)
    if(is.null(main.panel))
      main.panel <- names(x)
    else {
      if(!is.expression(main.panel))
        main.panel <- as.character(main.panel)
      nmp <- length(main.panel)
      if(nmp == 1)
        main.panel <- rep.int(main.panel, n)
      else if(nmp != n)
        stop("Incorrect length for main.panel")
    }

    if(allfv && equal.scales) {
      ## all entries are 'fv' objects: determine their plot limits
      fvlims <- lapply(x, plot, ..., limitsonly=TRUE)
      ## establish common x,y limits for all panels
      xlim <- range(unlist(lapply(fvlims, getElement, name="xlim")))
      ylim <- range(unlist(lapply(fvlims, getElement, name="ylim")))
      extrargs <- list(xlim=xlim, ylim=ylim)
    } else extrargs <- list()
    
    if(!arrange) {
      ## sequence of plots
      for(i in 1:n) {
        xi <- x[[i]]
        exec.or.plot(panel.begin, i, xi, main=main.panel[i],
                     extrargs=extrargs)
        extraplot(i, xi, ...,
                  add=!is.null(panel.begin),
                  main=main.panel[i],
                  panel.args=panel.args, extrargs=extrargs,
                  plotcommand=plotcommand)
        exec.or.plot(panel.end, i, xi, add=TRUE, extrargs=extrargs)
      }
      if(!is.null(adorn.left))
        warning("adorn.left was ignored because arrange=FALSE")
      if(!is.null(adorn.right))
        warning("adorn.right was ignored because arrange=FALSE")
      if(!is.null(adorn.top))
        warning("adorn.top was ignored because arrange=FALSE")
      if(!is.null(adorn.bottom))
        warning("adorn.bottom was ignored because arrange=FALSE")
      return(invisible(NULL))
    }

    ## ARRAY of plots
    ## decide whether to plot a main header
    main <- if(!missing(main) && !is.null(main)) main else xname
    if(!is.character(main)) {
      ## main title could be an expression
      nlines <- 1
      banner <- TRUE
    } else {
      ## main title is character string/vector, possibly ""
      banner <- any(nzchar(main))
      if(length(main) > 1)
        main <- paste(main, collapse="\n")
      nlines <- length(unlist(strsplit(main, "\n")))
    }
    ## determine arrangement of plots
    ## arrange like mfrow(nrows, ncols) plus a banner at the top
    if(is.null(nrows) && is.null(ncols)) {
      nrows <- as.integer(floor(sqrt(n)))
      ncols <- as.integer(ceiling(n/nrows))
    } else if(!is.null(nrows) && is.null(ncols))
      ncols <- as.integer(ceiling(n/nrows))
    else if(is.null(nrows) && !is.null(ncols))
      nrows <- as.integer(ceiling(n/ncols))
    else stopifnot(nrows * ncols >= length(x))
    nblank <- ncols * nrows - n
    if(allfv || list(plotcommand) %in% list("persp", persp)) {
      ## Function plots do not have physical 'size'
      sizes.known <- FALSE
    } else {
      ## Determine dimensions of objects
      ##     (including space for colour ribbons, if they are images)
      boxes <- getPlotBoxes(x, ..., plotcommand=plotcommand,
                            panel.args=panel.args, extrargs=extrargs)
      sizes.known <- !any(sapply(boxes, inherits, what="try-error"))
      if(sizes.known) {
        extrargs <- resolve.defaults(extrargs, list(claim.title.space=TRUE))
        boxes <- getPlotBoxes(x, ..., plotcommand=plotcommand,
                              panel.args=panel.args, extrargs=extrargs)
      }
      if(equal.scales && !sizes.known) {
        warning("Ignored equal.scales=TRUE; scales could not be determined")
        equal.scales <- FALSE
      }
    }
    if(sizes.known) {
      ## determine size of each panel
      if(equal.scales) {
        ## do not rescale panels
        scaledboxes <- boxes
      } else {
        ## rescale panels
        sides <- lapply(boxes, sidelengths)
        bwidths <- unlist(lapply(sides, "[", 1))
        bheights <- unlist(lapply(sides, "[", 2))
        ## Force equal heights, unless there is only one column
        scales <- if(ncols > 1) 1/bheights else 1/bwidths
        scaledboxes <- vector(mode="list", length=n)
        for(i in 1:n)
          scaledboxes[[i]] <- scalardilate(boxes[[i]], scales[i])
      }
    }
    ## determine whether to display all objects in one enormous plot
    ## Precondition is that everything has a spatial bounding box
    single.plot <- equal.scales && sizes.known
    if(equal.scales && !single.plot && !allfv)
      warning("equal.scales=TRUE ignored ", "because bounding boxes ",
              "could not be determined", call.=FALSE)
    ## enforce alignment by expanding boxes
    if(halign) {
      if(!equal.scales)
        warning("halign=TRUE ignored because equal.scales=FALSE")
      ## x coordinates align in each column
      xr <- range(sapply(scaledboxes, getElement, name="xrange"))
      scaledboxes <- lapply(scaledboxes, "[[<-", i="xrange", value=xr)
    }
    if(valign) {
      if(!equal.scales)
        warning("valign=TRUE ignored because equal.scales=FALSE")
      ## y coordinates align in each column
      yr <- range(sapply(scaledboxes, getElement, name="yrange"))
      scaledboxes <- lapply(scaledboxes, "[[<-", i="yrange", value=yr)
    }
    ## set up layout
    mat <- matrix(c(seq_len(n), integer(nblank)),
                  byrow=TRUE, ncol=ncols, nrow=nrows)
    if(sizes.known) {
      boxsides <- lapply(scaledboxes, sidelengths)
      xwidths <- unlist(lapply(boxsides, "[", i=1))
      xheights <- unlist(lapply(boxsides, "[", i=2))
      heights <- apply(mat, 1, function(j,h) { max(h[j[j>0]]) }, h=xheights)
      widths <- apply(mat, 2, function(i,w) { max(w[i[i>0]]) }, w=xwidths)
    } else {
      heights <- rep.int(1, nrows)
      widths <- rep.int(1, ncols)
    }
    meanheight <- mean(heights)
    meanwidth  <- mean(widths)
    nall <- n
    ##
    if(single.plot) {
      ## .........  create a single plot ..................
      ## determine sizes
      ht <- max(heights)
      wd <- max(widths)
      marpar <- mar.panel * c(ht, wd, ht, wd)/6
      vsep <- vsep * ht/6
      hsep <- hsep * wd/6
      mainheight <- any(nzchar(main.panel)) * ht/5
      ewidths <- marpar[2] + widths + marpar[4]
      eheights <- marpar[1] + heights + marpar[3] + mainheight
      Width <- sum(ewidths) + hsep * (length(ewidths) - 1)
      Height <- sum(eheights) + vsep * (length(eheights) - 1)
      bigbox <- owin(c(0, Width), c(0, Height))
      ox <- marpar[2] + cumsum(c(0, ewidths + hsep))[1:ncols]
      oy <- marpar[1] + cumsum(c(0, rev(eheights) + vsep))[nrows:1]
      panelorigin <- as.matrix(expand.grid(x=ox, y=oy))
      ## initialise, with banner
      cex <- resolve.1.default(list(cex.title=1.5), list(...))/par('cex.main')
      plot(bigbox, type="n", main=main, cex.main=cex)
      ## plot individual objects
      for(i in 1:n) {
        ## determine shift vector that moves bottom left corner of spatial box
        ## to bottom left corner of target area on plot device
        vec <- panelorigin[i,] - with(scaledboxes[[i]], c(xrange[1], yrange[1]))
        ## shift panel contents
        xi <- x[[i]]
        xishift <- shift(xi, vec)
        ## let rip
        if(!is.null(panel.begin))
          exec.or.plotshift(panel.begin, i, xishift,
                            add=TRUE,
                            main=main.panel[i], show.all=TRUE,
                            extrargs=extrargs,
                            vec=vec)
        extraplot(i, xishift, ...,
                  add=TRUE, show.all=is.null(panel.begin),
                  main=main.panel[i],
                  extrargs=extrargs,
                  panel.args=panel.args, plotcommand=plotcommand)
        exec.or.plotshift(panel.end, i, xishift, add=TRUE,
                          extrargs=extrargs,
                          vec=vec)
      }
      return(invisible(NULL))
    }
    ## ................. multiple logical plots using 'layout' ..............
    ## adjust panel margins to accommodate desired extra separation
    mar.panel <- pmax(0, mar.panel + c(vsep, hsep, vsep, hsep)/2)
    ## check for adornment
    if(!is.null(adorn.left)) {
      ## add margin at left, of width adorn.size * meanwidth
      nall <- i.left <- n+1
      mat <- cbind(i.left, mat)
      widths <- c(adorn.size * meanwidth, widths)
    } 
    if(!is.null(adorn.right)) {
      ## add margin at right, of width adorn.size * meanwidth
      nall <- i.right <- nall+1
      mat <- cbind(mat, i.right)
      widths <- c(widths, adorn.size * meanwidth)
    } 
    if(!is.null(adorn.bottom)) {
      ## add margin at bottom, of height adorn.size * meanheight
      nall <- i.bottom <- nall+1
      mat <- rbind(mat, i.bottom)
      heights <- c(heights, adorn.size * meanheight)
    } 
    if(!is.null(adorn.top)) {
      ## add margin at top, of height adorn.size * meanheight
      nall <- i.top <- nall + 1
      mat <- rbind(i.top, mat)
      heights <- c(adorn.size * meanheight, heights)
    } 
    if(banner) {
      ## Increment existing panel numbers
      ## New panel 1 is the banner
      panels <- (mat > 0)
      mat[panels] <- mat[panels] + 1
      mat <- rbind(1, mat)
      heights <- c(0.1 * meanheight * (1 + nlines), heights)
    }
    ## declare layout
    layout(mat, heights=heights, widths=widths, respect=sizes.known)
    ## start output .....
    ## .... plot banner
    if(banner) {
      opa <- par(mar=rep.int(0,4), xpd=TRUE)
      plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
           xlim=c(-1,1),ylim=c(-1,1))
      cex <- resolve.1.default(list(cex.title=1.5), list(...))/par('cex')
      text(0,0,main, cex=cex)
    }
    ## plot panels
    npa <- par(mar=mar.panel)
    if(!banner) opa <- npa
    for(i in 1:n) {
      xi <- x[[i]]
      exec.or.plot(panel.begin, i, xi, main=main.panel[i], extrargs=extrargs)
      extraplot(i, xi, ...,
                add = !is.null(panel.begin), 
                main = main.panel[i],
                extrargs=extrargs,
                panel.args=panel.args, plotcommand=plotcommand)
      exec.or.plot(panel.end, i, xi, add=TRUE, extrargs=extrargs)
    }
    ## adornments
    if(nall > n) {
      par(mar=rep.int(0,4), xpd=TRUE)
      if(!is.null(adorn.left))
        adorn.left()
      if(!is.null(adorn.right))
        adorn.right()
      if(!is.null(adorn.bottom))
        adorn.bottom()
      if(!is.null(adorn.top))
        adorn.top()
    }
    ## revert
    layout(1)
    par(opa)
    return(invisible(NULL))
  }

  plot.anylist
})


contour.imlist <- contour.listof <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  do.call("plot.solist",
          resolve.defaults(list(x=x, plotcommand="contour"),
                           list(...),
                           list(main=xname)))
}

plot.imlist <- local({

  plot.imlist <- function(x, ..., plotcommand="image",
                          equal.ribbon = FALSE, ribmar=NULL) {
    xname <- short.deparse(substitute(x))
    if(equal.ribbon &&
       (list(plotcommand) %in% list("image", "plot", image, plot))) {
      out <- imagecommon(x, ..., xname=xname, ribmar=ribmar)
    } else {
      out <- do.call("plot.solist",
                     resolve.defaults(list(x=x, plotcommand=plotcommand), 
                                      list(...),
                                      list(main=xname)))
    }
    return(invisible(out))
  }

  imagecommon <- function(x, ...,
                          xname,
                          zlim=NULL,
                          ribbon=TRUE,
                          ribside=c("right", "left", "bottom", "top"),
                          ribsep=NULL, ribwid=0.5, ribn=1024,
                          ribscale=NULL, ribargs=list(),
                          ribmar = NULL, mar.panel = c(2,1,1,2)) {
    if(missing(xname))
      xname <- short.deparse(substitute(x))
    ribside <- match.arg(ribside)
    stopifnot(is.list(ribargs))
    if(!is.null(ribsep))
      warning("Argument ribsep is not yet implemented for image arrays")
    ## determine range of values
    if(is.null(zlim))
      zlim <- range(unlist(lapply(x, range)))
    ## determine common colour map
    imcolmap <- plot.im(x[[1]], do.plot=FALSE, zlim=zlim, ..., ribn=ribn)
    ## plot ribbon?
    if(!ribbon) {
      ribadorn <- list()
    } else {
      ## determine plot arguments for colour ribbon
      vertical <- (ribside %in% c("right", "left"))
      scaleinfo <- if(!is.null(ribscale)) list(labelmap=ribscale) else list()
      sidecode <- match(ribside, c("bottom", "left", "top", "right"))
      ribstuff <- c(list(x=imcolmap, main="", vertical=vertical),
                    ribargs,
                    scaleinfo,
                    list(side=sidecode))
      if (is.null(mar.panel)) 
        mar.panel <- c(2, 1, 1, 2)
      if (length(mar.panel) != 4) 
        mar.panel <- rep(mar.panel, 4)[1:4]
      if (is.null(ribmar)) {
        ribmar <- mar.panel/2
        newmar <- c(2, 0)
        switch(ribside,
               left   = { ribmar[c(2, 4)] <- newmar },
               right  = { ribmar[c(4, 2)] <- newmar },
               bottom = { ribmar[c(1, 3)] <- newmar },
               top    = { ribmar[c(3, 1)] <- newmar }
               )
      }
       ## function executed to plot colour ribbon
      do.ribbon <- function() {
        opa <- par(mar=ribmar)
        do.call("plot", ribstuff)
        par(opa)
      }
      ## encoded as 'adorn' argument
      ribadorn <- list(adorn=do.ribbon, adorn.size=ribwid)
      names(ribadorn)[1] <- paste("adorn", ribside, sep=".")
    }
    ##
    do.call("plot.solist",
            resolve.defaults(list(x=x, plotcommand="image"),
                             list(...),
                             list(mar.panel=mar.panel,
                                  main=xname,
                                  col=imcolmap, zlim=zlim,
                                  ribbon=FALSE),
                             ribadorn))
  }

  plot.imlist
})

image.imlist <- image.listof <-
  function(x, ..., equal.ribbon = FALSE, ribmar=NULL) {
    plc <- resolve.1.default(list(plotcommand="image"), list(...))
    if(list(plc) %in% list("image", "plot", image, plot)) {
      out <- plot.imlist(x, ..., plotcommand="image",
                         equal.ribbon=equal.ribbon, ribmar=ribmar)
    } else {
      out <- plot.solist(x, ..., ribmar=ribmar)
    }
    return(invisible(out))
  }

